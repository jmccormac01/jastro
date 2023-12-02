"""
Functions for reducing data
"""
import numpy as np
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
        master_bias, _ = jhk.load_fits_image(master_bias_filename)
        return master_bias
    except:
        # check for no images
        if not images.files:
            return None
        for f in images.files_filtered(imagetyp=bias_keyword):
            print(f)
            ccd, _ = jhk.load_fits_image(f)
            bias_list.append(ccd)
        bias_list = np.array(bias_list)

        if len(bias_list) > 0:
            new_header = {}
            new_header["N_STACK"] = len(bias_list)
            master_bias = np.median(bias_list, axis=0)
            jhk.write_fits_image(master_bias_filename, master_bias,
                                 header=new_header, clobber=True)
            return master_bias
        return None

def make_master_dark(images, master_bias=None, dark_keyword='DARK',
                     exptime_keyword='EXPTIME', master_dark_filename="master_dark.fits",
                     med_bias=0):
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
    med_bias : int
        Fallback value if overscan is unavailable
        Default = 0

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
        master_dark, hdr = jhk.load_fits_image(master_dark_filename)
        dark_exp = round(float(hdr[exptime_keyword]), 2)
        return master_dark, dark_exp
    except:
        # check for no images
        if not images.files:
            return None, None
        for f in images.files_filtered(imagetyp=dark_keyword):
            print(f)
            ccd, hdr = jhk.load_fits_image(f)
            if not dark_exp:
                dark_exp = round(float(hdr[exptime_keyword]), 2)
            if master_bias is not None:
                ccd = ccd - master_bias
            else:
                print(f'No master bias, subtracting med_bias {med_bias}...')
                ccd = ccd - med_bias
            dark_list.append(ccd)
        dark_list = np.array(dark_list)

        if len(dark_list) > 0:
            new_header = {}
            master_dark = np.median(dark_list, axis=0)
            new_header["N_STACK"] = len(dark_list)
            new_header[exptime_keyword] = dark_exp
            jhk.write_fits_image(master_dark_filename, master_dark,
                                 header=new_header, clobber=True)
            return master_dark, dark_exp
        return None, None

def make_master_dark_osc(images, overscan_keyword, dark_keyword='DARK',
                         exptime_keyword='EXPTIME', master_dark_filename="master_dark.fits",
                         med_bias=0):
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
    med_bias : int
        Fallback value if overscan is unavailable
        Default = 0

    Returns
    -------
    master_dark : array-like | None
        A master dark image array, or
        None if no darks are found
    dark_exp : float | None
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
        dark_exp = round(float(header[exptime_keyword]), 2)
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
            dark_exp = round(float(header[exptime_keyword]), 2)
            # fetch overscan correction
            try:
                os_region = header[overscan_keyword]
                os_corr = extract_overscan_correction(ccd, os_region)
            except KeyError:
                print(f'No overscan, subtracting med_bias {med_bias}...')
                os_corr = med_bias
            # correct the frame
            ccd_corr = ccd - os_corr
            # save dark signal for combining
            dark_list.append(ccd_corr)
        dark_list = np.array(dark_list)

        if len(dark_list) > 0:
            new_header = {}
            master_dark = np.median(dark_list, axis=0)
            new_header["N_STACK"] = len(dark_list)
            new_header[exptime_keyword] = dark_exp
            jhk.write_fits_image(master_dark_filename, master_dark,
                                 header=new_header, clobber=True)
            return master_dark, dark_exp
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
                     master_flat_filename="master_flat.fits", med_bias=0):
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
    med_bias : int
        Fallback value if overscan is unavailable
        Default = 0

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
        master_flat, _ = jhk.load_fits_image(master_flat_filename)
        return master_flat
    except:
        # check for no images
        if not images.files:
            return None
        # create the master flat field for each filter
        print(f'Reducing flats from filter {filt}')
        for f in images.files_filtered(imagetyp=flat_keyword, filter=filt):
            print(f)
            ccd, hdr = jhk.load_fits_image(f)
            data_exp = round(float(hdr[exptime_keyword]), 2)
            if master_bias is not None:
                ccd = ccd - master_bias
            else:
                print(f'No master bias, subtracting med_bias {med_bias}...')
                ccd = ccd - med_bias
            if master_dark is not None:
                ccd = ccd - (master_dark * (data_exp/dark_exp))
            else:
                print('No master dark, skipping correction...')
            sky_level, _ = estimate_sky_level(ccd)
            ccd = ccd / sky_level
            flat_list.append(ccd)
        flat_list = np.array(flat_list)

        if len(flat_list) > 0:
            master_flat = np.median(flat_list, axis=0)
            jhk.write_fits_image(master_flat_filename, master_flat,
                                 header=False, clobber=True)
            return master_flat
        print(f'There are no flats for {filt}, skipping...')
        return None

def make_master_flat_osc(images, filt, overscan_keyword, master_dark=None,
                         dark_exp=30, flat_keyword='FLAT', exptime_keyword='EXPTIME',
                         master_flat_filename="master_flat.fits", med_bias=0):
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
    overscan_keyword : str
        Overscan region defined as '[x1:x2,y1:y2]'
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
    med_bias : int
        Fallback value if overscan is unavailable
        Default = 0

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
        master_flat, _ = jhk.load_fits_image(master_flat_filename)
        return master_flat
    except:
        # check for no images
        if not images.files:
            return None
        # create the master flat field for each filter
        print(f'Reducing flats from filter {filt}')
        for f in images.files_filtered(imagetyp=flat_keyword, filter=filt):
            print(f)
            # load the data
            ccd, header = jhk.load_fits_image(f)
            # load few header items
            data_exp = round(float(header[exptime_keyword]), 2)
            # fetch overscan correction
            try:
                os_region = header[overscan_keyword]
                os_corr = extract_overscan_correction(ccd, os_region)
            except KeyError:
                print(f'No master bias, subtracting med_bias {med_bias}...')
                os_corr = med_bias
            # correct the frame for the overscan
            ccd_corr = ccd - os_corr
            # correct for dark current if master_dark
            if master_dark is not None:
                # scale the dark signal before correcting
                ccd_corr = ccd_corr - (master_dark * (data_exp/dark_exp))
            else:
                print('No master dark, skipping correction...')
            # estimate sky level
            sky_level, _ = estimate_sky_level(ccd_corr)
            # normalise by the sky level
            ccd_corr_norm = ccd_corr/sky_level
            # save dark signal for combining
            flat_list.append(ccd_corr_norm)
        flat_list = np.array(flat_list)

        if len(flat_list) > 0:
            master_flat = np.median(flat_list, axis=0)
            jhk.write_fits_image(master_flat_filename, master_flat,
                                 header=False, clobber=True)
            return master_flat
        print(f'There are no flats for {filt}, skipping...')
        return None

def correct_data(filename, filt, location, master_bias=None, master_dark=None,
                 master_flat=None, dark_exp=30, exptime_keyword='EXPTIME',
                 dateobs_start_keyword='DATE-OBS', ra_keyword='RA', dec_keyword='DEC',
                 med_bias=0):
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
    med_bias : int
        Fallback value if overscan is unavailable
        Default = 0

    Returns
    -------
    ccd :  array
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
    ccd, hdr = jhk.load_fits_image(filename)
    # collect/correct some header values
    data_exp = round(float(hdr[exptime_keyword]), 2)
    half_exptime = data_exp/2.
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

    if master_bias is not None:
        ccd = ccd - master_bias
    else:
        print(f'No master bias, subtracting med_bias {med_bias}...')
        ccd = ccd - med_bias
    if master_dark is not None:
        ccd = ccd - (master_dark * (data_exp/dark_exp))
    else:
        print('No master dark, skipping correction...')
    if master_flat is not None:
        ccd = ccd / master_flat
    else:
        print(f'No master flat for {filt}, skipping correction...')

    # after calibrating we need np.float64 data
    # if there are no calibrations we maintain dtype = np.uint16
    # fix this by doing the following
    if isinstance(ccd[0][0], np.uint16):
        ccd = ccd.astype(np.float64)
    return ccd, time_jd, time_bary, time_helio

def correct_data_osc(filename, filt, location, overscan_keyword, master_dark=None,
                     master_flat=None, dark_exp=30, exptime_keyword='EXPTIME',
                     dateobs_start_keyword='DATE-OBS', ra_keyword='RA', dec_keyword='DEC',
                     med_bias=0):
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
    overscan_keyword : str
        Overscan region defined as '[x1:x2,y1:y2]'
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
    med_bias : int
        Fallback value if overscan is unavailable
        Default = 0

    Returns
    -------
    ccd : array
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

    # load the data
    ccd, header = jhk.load_fits_image(filename)
    # load few header items
    data_exp = round(float(header[exptime_keyword]), 2)
    ra = header[ra_keyword]
    dec = header[dec_keyword]
    # do some time conversion
    half_exptime = data_exp/2.
    time_isot = Time(header[dateobs_start_keyword], format='isot',
                     scale='utc', location=location)
    time_jd = Time(time_isot.jd, format='jd', scale='utc', location=location)
    # correct to mid exposure time
    time_jd = time_jd + half_exptime*u.second
    ltt_bary, ltt_helio = jcoords.get_light_travel_times(ra, dec, time_jd)
    time_bary = time_jd.tdb + ltt_bary
    time_helio = time_jd.utc + ltt_helio

    # fetch overscan correction
    try:
        os_region = header[overscan_keyword]
        os_corr = extract_overscan_correction(ccd, os_region)
    except KeyError:
        print(f'No master bias, subtracting med_bias {med_bias}...')
        # make a dummy block of zeros for OSC
        os_corr = med_bias

    # correct the frame for the overscan
    ccd_corr = ccd - os_corr
    # correct for dark current if master_dark
    if master_dark is not None:
        # scale the dark signal before correcting
        ccd_corr = ccd_corr - (master_dark * (data_exp/dark_exp))
    else:
        print('No master dark, skipping correction...')
    if master_flat is not None:
        ccd_corr = ccd_corr / master_flat
    else:
        print(f'No master flat for {filt}, skipping correction...')
    # check the data type is float for sep stage
    if isinstance(ccd_corr[0][0], np.uint16):
        ccd_corr = ccd_corr.astype(np.float64)
    return ccd_corr, time_jd, time_bary, time_helio

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
