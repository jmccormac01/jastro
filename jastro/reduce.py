"""
Functions for reducing data
"""
from collections import defaultdict
import numpy as np
import ccdproc
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import jastro.coords as jcoords

# pylint: disable=invalid-name
# pylint: disable=no-member

def get_image_list(directory='.', glob_exclude='master*'):
    """
    Make a list of what images are present in a given
    directory

    Parameters
    ----------
    directory : str
        The path to the folder to analyse for images
        Default = '.' (the current directory)

    Returns
    -------
        The ccdproc ImageFileCollection of that directory

    Raises
    ------
    None
    """
    return ccdproc.ImageFileCollection(directory, glob_exclude=glob_exclude)

def make_master_bias(images, bias_keyword='BIAS'):
    """
    If no master bias image is found try making a
    master bias using all biases found in the
    ImageFileCollection object, if any.

    Parameters
    ----------
    images : object ccdproc.ImageFileCollection
        Object containing a list of images in the
        working directory
    bias_keyword : str
        Header keyword for bias images
        Default = 'BIAS'

    Returns
    -------
    master_bias : array-like | None
        A master bias image array, or
        None if no biases are found

    Raises
    ------
    None
    """
    bias_list = []
    try:
        master_bias = ccdproc.CCDData.read('master_bias.fits', unit=u.adu)
        return master_bias
    except FileNotFoundError:
        # check for no images
        if not images.files:
            return None
        for f in images.files_filtered(imagetyp=bias_keyword):
            print(f)
            ccd = ccdproc.CCDData.read(f, unit=u.adu)
            bias_list.append(ccd)
    try:
        master_bias = ccdproc.combine(bias_list, method='median')
        master_bias.write('master_bias.fits', overwrite=True)
        return master_bias
    except IndexError:
        return  None

def make_master_dark(images, master_bias=None, dark_keyword='DARK',
                     exptime_keyword='EXPTIME'):
    """
    If no master dark image is found try making a
    master dark from all darks found in the
    ImageFileCollection object.

    If a master bias image is provided the darks
    are first corrected for their bias level
    before combination

    Parameters
    ----------
    images : object ccdproc.ImageFileCollection
        Object containing a list of images in the
        working directory
    master_bias : array-like, optional
        Array containing the master bias image
        Default = None
    dark_keyword : str
        Header keyword for dark images
        Default = 'DARK'
    exptime_keyword : str
        Header keyword for exposure time
        Default = 'EXPTIME'

    Returns
    -------
    master_dark : array-like | None
        A master dark image array, or
        None if no darks are found
    dark_exp : int | None
        The exposure time of the dark frames, or
        None if no darks are found

    Raises
    ------
    None
    """
    dark_list = []
    dark_exp = None
    try:
        fitsfile = 'master_dark.fits'
        master_dark = ccdproc.CCDData.read(fitsfile, unit=u.adu)
        dark_exp = int(fits.open(fitsfile)[0].header[exptime_keyword])
        return master_dark, dark_exp
    except FileNotFoundError:
        # check for no images
        if not images.files:
            return None, None
        for f in images.files_filtered(imagetyp=dark_keyword):
            print(f)
            if not dark_exp:
                with fits.open(f) as fitsfile:
                    dark_exp = int(fitsfile[0].header[exptime_keyword])
            ccd = ccdproc.CCDData.read(f, unit=u.adu)
            if master_bias:
                ccd = ccdproc.subtract_bias(ccd, master_bias)
            else:
                print('No master bias, skipping correction...')
            dark_list.append(ccd)
    try:
        master_dark = ccdproc.combine(dark_list, method='median')
        master_dark.write('master_dark.fits', overwrite=True)
        return master_dark, dark_exp
    except IndexError:
        return None, None

def estimate_sky_level(data):
    """
    Function to interatively sigma clip the sky background
    to estimate the sky level without the influence of stars
    """
    mean_diff = 1E6
    mean_diff_limit = 1E-6
    sigma = 3
    # create a masked array where nothing is masked
    data = np.ma.masked_where(data < -1E6, data)
    while mean_diff > mean_diff_limit:
        mean = np.ma.average(data)
        rms = np.ma.std(data)
        masked_data = np.ma.masked_where(((data > mean+sigma*rms) | (data < mean-sigma*rms)), data)
        new_mean = np.ma.average(masked_data)
        new_rms = np.ma.std(masked_data)
        print('Sky level: {}, RMS: {}'.format(new_mean, new_rms))
        data = masked_data
        mean_diff = abs(new_mean-mean)/new_mean
    return new_mean, new_rms

def make_master_flat(images, filt, master_bias=None, master_dark=None,
                     dark_exp=30, flat_keyword='Flat Field',
                     exptime_keyword='EXPTIME'):
    """
    If no master flat is found try making a master flat
    from all the flats in the ImageFileCollection
    object.

    If a master bias is provided the flats are first
    corrected for their bias level.

    Similarly, if a master dark is provided the flats
    are corrected for their dark current.

    Parameters
    ----------
    images : object ccdproc.ImageFileCollection
        Object containing a list of images in the
        working directory
    filt : str
        The name of the filter used for this data
    master_bias : array-like, optional
        Array containing the master bias image
        Default = None
    master_dark : array-like, optional
        Array containing the master dark image
        Default = None
    dark_exp : int, optional
        Exposure time for master dark
        Default = 30
    flat_keyword : str, optional
        Header keyword for the flat fields
        Default = 'Flat Field'
    exptime_keyword : str, optional
        Header keyword for exposure time
        Default = 'EXPTIME'

    Returns
    -------
    master_flat : array-like | None
        A master flat image array, or
        None if no flats are found

    Raises
    ------
    None
    """
    # empty dictionaries for the filtered data
    flat_list = defaultdict(list)
    try:
        fitsfile = 'master_flat_{0:s}.fits'.format(filt)
        master_flat = ccdproc.CCDData.read(fitsfile, unit=u.adu)
        return master_flat
    except FileNotFoundError:
        # check for no images
        if not images.files:
            return None
        # create the master flat field for each filter
        print('Reducing flats from filter {0:s}'.format(filt))
        for f in images.files_filtered(imagetyp=flat_keyword, filter=filt):
            print(f)
            with fits.open(f) as fitsfile:
                data_exp = int(fitsfile[0].header[exptime_keyword])
            ccd = ccdproc.CCDData.read(f, unit=u.adu)
            if master_bias:
                ccd = ccdproc.subtract_bias(ccd, master_bias)
            else:
                print('No master bias, skipping correction...')
            if master_dark:
                ccd = ccdproc.subtract_dark(ccd, master_dark,
                                            scale=True,
                                            dark_exposure=dark_exp*u.second,
                                            data_exposure=data_exp*u.second)
            else:
                print('No master dark, skipping correction...')
            sky_level, _ = estimate_sky_level(ccd.data)
            ccd.data = ccd.data/sky_level
            flat_list[filt].append(ccd)
    try:
        master_flat = ccdproc.combine(flat_list[filt], method='median')
        master_flat.write('master_flat_{0:s}.fits'.format(filt), overwrite=True)
        return master_flat
    except IndexError:
        print('There are no flats for {0:s}, skipping...'.format(filt))
        master_flat = None

def correct_data(filename, filt, location, master_bias=None, master_dark=None,
                 master_flat=None, dark_exp=30, exptime_keyword='EXPTIME',
                 jd_keyword='JD', ra_keyword='RA', dec_keyword='DEC',
                 output_reduced_frames=False):
    """
    Correct a science image using the available
    master calibrations. Skip a calibration step if the
    master frame does not exist.

    No reduced file is written in this new scheme.
    Instead, the corrected data is passed directly
    to the phot() routine, photometry is done as per
    the configuration and the photometry is written out
    only.

    Parameters
    ----------
    filename : str
        Name of the image to process
    filt : str
        The name of the filter used for this data
    location : astropy.coordinates.EarthLocation
        The EarthLocation of the observatory where the data
        was taken. This is used for light travel time
        calculations
    master_bias : array-like, optional
        Array containing the master bias image
        Default = None
    master_dark : array-like, optional
        Array containing the master dark image
        Default = None
    master_flat : array-like, optional
        Array containing the master flat image
        Default = None
    dark_exp : int, optional
        Exposure time for master dark
        Default = 30
    exptime_keyword : str, optional
        Header keyword for exposure time
        Default = 'EXPTIME'
    jd_keyword : str, optional
        Header keyword for the Julian Date
        Default = 'JD'
    ra_keyword : str, optional
        Header keyword for the Right Ascension
        Default = 'RA'
    dec_keyword : str. optional
        Header keyword for the Declination
        Default = 'DEC'
    output_reduced_frames : bool, optional
        Output the reduced frames?
        Default = False

    Returns
    -------
    ccd : ccdproc.CCDData array
        The corrected image
    time_jd : astropy.Time
        The JD at mid exposure
    time_bary : astropy.Time
        The BJD_TDB at mid exposure
    time_helio : astropy.Time
        The HJD at mid exposure

    Raises
    ------
    None
    """
    print('Reducing {0:s}...'.format(filename))
    with fits.open(filename) as fitsfile:
        # collect/correct some header values
        hdr = fitsfile[0].header
        data_exp = int(hdr[exptime_keyword])
        half_exptime = float(data_exp)/2.
        time_jd = Time(float(hdr[jd_keyword]), format='jd',
                       scale='utc', location=location)
        # correct to mid exposure time
        time_jd = time_jd + half_exptime*u.second
        ra = hdr[ra_keyword]
        dec = hdr[dec_keyword]
        ltt_bary, ltt_helio = jcoords.get_light_travel_times(ra, dec, time_jd)
        time_bary = time_jd.tdb + ltt_bary
        time_helio = time_jd.utc + ltt_helio

    ccd = ccdproc.CCDData.read(filename, unit=u.adu)
    if master_bias:
        ccd = ccdproc.subtract_bias(ccd, master_bias)
    else:
        print('No master bias, skipping correction...')
    if master_dark:
        ccd = ccdproc.subtract_dark(ccd, master_dark,
                                    scale=True,
                                    dark_exposure=dark_exp*u.second,
                                    data_exposure=data_exp*u.second)
    else:
        print('No master dark, skipping correction...')
    if master_flat:
        ccd = ccdproc.flat_correct(ccd, master_flat)
    else:
        print('No master flat for {0:s}, skipping correction...'.format(filt))

    # after calibrating we get np.float64 data
    # if there are no calibrations we maintain dtype = np.uint16
    # sep weeps
    # fix this by doing the following
    if isinstance(ccd.data[0][0], np.uint16):
        ccd.data = ccd.data.astype(np.float64)
    # output the files
    if output_reduced_frames:
        new_filename = '{}_r.fts'.format(filename.split('.fts')[0])
        fits.writeto(new_filename, ccd.data, header=hdr, overwrite=True)
    return ccd, time_jd, time_bary, time_helio
