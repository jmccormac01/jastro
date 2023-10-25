"""
Functions for reducing data
"""
import numpy as np
import ccdproc
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import jastro.coords as jcoords
import jastro.housekeeping as jhk

# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=bare-except

def extract_overscan_correction(data, os_region):
    """
    Take an image and extract the overscan correction
    """
    x1, x2 = map(int, os_region.split(',')[0][1:].split(":"))
    y1, y2 = map(int, os_region.split(',')[1][:-1].split(":"))

    # extract the overscan region
    os_data = data[y1-1:y2, x1-1:x2]

    # detemine direction of overscan, along row or columns?
    if x2-x1 > y2-y1:
        os_correction = np.median(os_data, axis=0)
    else:
        os_correction = np.median(os_data, axis=1)
    return os_correction

def make_master_bias(images, bias_keyword='BIAS',
                     master_bias_filename="master_bias.fits"):
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
    master_bias_filename : str
        Name of master bias file to save
        Default = 'master_bias.fits'

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
        master_bias = ccdproc.CCDData.read(master_bias_filename, unit=u.adu)
        return master_bias
    except:
        # check for no images
        if not images.files:
            return None
        for f in images.files_filtered(imagetyp=bias_keyword):
            print(f)
            ccd = ccdproc.CCDData.read(f, unit=u.adu)
            bias_list.append(ccd)
    try:
        master_bias = ccdproc.combine(bias_list, method='median')
        master_bias.write(master_bias_filename, overwrite=True)
        return master_bias
    except IndexError:
        return  None

def make_master_dark(images, master_bias=None, dark_keyword='DARK',
                     exptime_keyword='EXPTIME', master_dark_filename="master_dark.fits"):
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
    master_dark_filename : str
        Name of master dark file to save
        Default = 'master_dark.fits'

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
        master_dark = ccdproc.CCDData.read(master_dark_filename, unit=u.adu)
        dark_exp = int(fits.open(master_dark_filename)[0].header[exptime_keyword])
        return master_dark, dark_exp
    except:
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
        master_dark.write(master_dark_filename, overwrite=True)
        return master_dark, dark_exp
    except IndexError:
        return None, None

def make_master_dark_osc(images, overscan_keyword, dark_keyword='DARK',
                         exptime_keyword='EXPTIME', master_dark_filename="master_dark.fits"):
    """
    If no master dark image is found try making a
    master dark from all darks found in the
    ImageFileCollection object.

    If a master bias image is provided the darks
    are first corrected for their bias level
    before combination

    Bias correction is done using the overscane region

    Parameters
    ----------
    images : object ccdproc.ImageFileCollection
        Object containing a list of images in the
        working directory
    overscan_keyword : str
        Overscan region defined as '[x1:x2,y1:y2]'
    dark_keyword : str
        Header keyword for dark images
        Default = 'DARK'
    exptime_keyword : str
        Header keyword for exposure time
        Default = 'EXPTIME'
    master_dark_filename : str
        Name of master dark file to save
        Default = 'master_dark.fits'

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
        master_dark, header = jhk.load_fits_image(master_dark_filename)
        dark_exp = int(header[exptime_keyword])
        return master_dark, dark_exp
    except:
        # check for no images
        if not images.files:
            return None, None
        for f in images.files_filtered(imagetyp=dark_keyword):
            print(f)
            # load the data
            ccd, header = jhk.load_fits_image(f)
            # load few header items
            dark_exp = int(header[exptime_keyword])
            os_region = header[overscan_keyword]
            # fetch overscan correction
            os_corr = extract_overscan_correction(ccd, os_region)
            # correct the frame
            ccd_corr = ccd - os_corr
            # save dark signal for combining
            dark_list.append(ccd_corr)
        dark_list = np.array(dark_list)

        if len(dark_list) > 0:
            master_dark = np.median(dark_list, axis=0)
            header["N_STACK"] = len(dark_list)
            jhk.write_fits_image(master_dark_filename, master_dark,
                                 header, clobber=True)
            return master_dark, dark_exp
        else:
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
        print(f'Sky level: {new_mean}, RMS: {new_rms}')
        data = masked_data
        mean_diff = abs(new_mean-mean)/new_mean
    return new_mean, new_rms

def make_master_flat(images, filt, master_bias=None, master_dark=None,
                     dark_exp=30, flat_keyword='FLAT', exptime_keyword='EXPTIME',
                     master_flat_filename="master_flat.fits"):
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
    master_flat_filename : str
        Name of master flat file to save
        Default = 'master_flat.fits'

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
    flat_list = []
    try:
        master_flat = ccdproc.CCDData.read(master_flat_filename, unit=u.adu)
        return master_flat
    except:
        # check for no images
        if not images.files:
            return None
        # create the master flat field for each filter
        print(f'Reducing flats from filter {filt}')
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
            flat_list.append(ccd)
    try:
        master_flat = ccdproc.combine(flat_list, method='median')
        master_flat.write(master_flat_filename, overwrite=True)
    except IndexError:
        print(f'There are no flats for {filt}, skipping...')
        master_flat = None
    return master_flat

def correct_data(filename, filt, location, master_bias=None, master_dark=None,
                 master_flat=None, dark_exp=30, exptime_keyword='EXPTIME',
                 dateobs_start_keyword='DATE-OBS', ra_keyword='RA', dec_keyword='DEC',
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
    dateobs_start_keyword : str, optional
        Header keyword for the date obs start
        Default = 'DATE-OBS'
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
    print(f'Reducing {filename}...')
    with fits.open(filename) as fitsfile:
        # collect/correct some header values
        hdr = fitsfile[0].header
        data_exp = int(hdr[exptime_keyword])
        half_exptime = float(data_exp)/2.
        time_isot = Time(hdr[dateobs_start_keyword], format='isot',
                         scale='utc', location=location)
        time_jd = Time(time_isot.jd, format='jd', scale='utc', location=location)
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
        print(f'No master flat for {filt}, skipping correction...')

    # after calibrating we get np.float64 data
    # if there are no calibrations we maintain dtype = np.uint16
    # sep weeps
    # fix this by doing the following
    if isinstance(ccd.data[0][0], np.uint16):
        ccd.data = ccd.data.astype(np.float64)
    # output the files
    if output_reduced_frames:
        prefix = filename.split('.fts')[0]
        new_filename = f'{prefix}_r.fts'
        fits.writeto(new_filename, ccd.data, header=hdr, overwrite=True)
    return ccd, time_jd, time_bary, time_helio

def find_max_pixel_value(data, x, y, radius):
    """
    Find the maximum pixel value in the image
    in a square around the aperture centre

    Parameters
    ----------
    data : array-like
        The image to search
    x : int
        X coordinate of the search box
    y : int
        Y coordinate of the search box
    radius : int
        The half width of the search box

    Returns
    -------
    max_pixel_value : int
        The maximum pixel value in the area provided

    Raises
    ------
    None
    """
    return round(data[int(y-radius):int(y+radius),
                      int(x-radius):int(x+radius)].ravel().max(), 2)
