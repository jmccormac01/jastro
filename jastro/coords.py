"""
Functions for handling on-sky or on chip coordinates
"""
import sys
import math
from datetime import datetime
import sep
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

# pylint: disable=invalid-name
# pylint: disable=no-member

def correct_proper_motion(night, pm_ra, pm_dec, coord_to_correct, epoch=2000):
    """
    Correct proper motion on a SkyCoord

    Parameters
    ----------
    night : str
        Night of the observation
    pm_ra : float
        mas/yr proper motion in RA
    pm_dec : float
        mas/yr proper motion in Dec
    coord_to_correct : astropy.coordinates.SkyCoord
        SkyCoord object to correct for proper motion

    Returns
    -------
    corrected_coord : astropy.coordinates.SkyCoord
        SkyCoord object corrected for proper motion

    Raises
    ------
    None
    """
    # find the number of decimal years since 2000.0
    # we ignore the fractional day as this is negligible
    night = datetime.strptime(night, "%Y%m%d").timetuple()
    year = night.tm_year
    year_day = night.tm_yday
    decimal_years_since_2000 = (year+(year_day/365.25)) - epoch
    # get the proper motion in degrees since 2000
    dec_cor = math.cos(math.radians(coord_to_correct.dec.deg))
    pm_ra_total = (pm_ra/(1000.0*60.*60.*dec_cor))*decimal_years_since_2000
    pm_dec_total = (pm_dec/(1000.0*60.*60.))*decimal_years_since_2000
    # correct the SkyCoord object
    corrected_coord = SkyCoord(ra=(coord_to_correct.ra.deg+pm_ra_total)*u.degree,
                               dec=(coord_to_correct.dec.deg+pm_dec_total)*u.degree,
                               frame='icrs')
    return corrected_coord

def catalogue_to_pixels(astrometry_image, catalogue_coords, object_ids):
    """
    Convert a list of catalogue positions to X and Y image
    coordinates

    Parameters
    ----------
    astrometry_image : str
        Name of the FITS file with solved WCS solution
    catalogue_coords : array-like
        RA and Dec in degrees of the targets positions to
        convert to pixels
    object_ids : array-like
        List of target matching IDs

    Returns
    -------
    x_checked : array-like
        X positions of stars found in the astrometry_image
    y_checked : array-like
        Y positions of stars found in the astrometry_image
    object_ids_checked : array-like
        Object IDs for the images found in the astrometry_image

    Raises
    ------
    None
    """
    try:
        with fits.open(astrometry_image) as fitsfile:
            hdr = fitsfile[0].header
            width = hdr['NAXIS1']
            height = hdr['NAXIS2']
            border = 50
    except FileNotFoundError:
        print('CANNOT FIND {}, EXITING...'.format(astrometry_image))
        sys.exit(1)

    # load the WCS
    w = WCS(hdr)
    # 0 is C indexing
    # 1 is Fortran indexing
    pix = w.wcs_world2pix(catalogue_coords, 0)
    x, y = pix[:, 0], pix[:, 1]
    # check that the pixel values are within the
    # bound of the image, exclude object if not
    x, y, object_ids = check_image_boundaries(x, y, object_ids, width, height, border)
    return x, y, object_ids

def pixels_to_catalogue(astrometry_image, pixel_coords):
    """
    Convert a list of X and Y image positions to coordinates

    Parameters
    ----------
    astrometry_image : str
        Name of the FITS file with solved WCS solution
    pixel_coords : array-like
        X and Y of the targets positions to
        convert to coordinates

    Returns
    -------
    ra : array-like
        RA positions in degrees of stars found in the astrometry_image
    dec : array-like
        DEC positions in degrees of stars found in the astrometry_image

    Raises
    ------
    None
    """
    try:
        with fits.open(astrometry_image) as fitsfile:
            hdr = fitsfile[0].header
    except FileNotFoundError:
        print('CANNOT FIND {}, EXITING...'.format(astrometry_image))
        sys.exit(1)

    # load the WCS
    w = WCS(hdr)
    # 0 is C indexing
    # 1 is Fortran indexing
    world = w.wcs_pix2world(pixel_coords, 0)
    ra, dec = world[:, 0], world[:, 1]
    return ra, dec

def check_image_boundaries(x, y, object_ids, width, height, border):
    """
    Check if a star is too close to the edge of an image

    Parameters
    ----------
    x : array-like
        X positions of stars to check
    y : array-like
        Y positions of stars to check
    object_ids : array-like
        List of object names for objects to check
    width : int
        CCD width in pixels
    height : int
        CCD height in pixels

    Returns
    -------
    x_checked : array-like
        X positions of objects safely on the image
    y_checked : array-like
        Y positions of objects safely on the image
    object_ids_checked : array-like
        Object IDs for objects safely on the image

    Raises
    ------
    None
    """
    x_checked, y_checked, object_ids_checked = [], [], []
    for i, j, k in zip(x, y, object_ids):
        if i > border and i < width-border:
            if j > border and j < height-border:
                x_checked.append(i)
                y_checked.append(j)
                object_ids_checked.append(k)
    return x_checked, y_checked, object_ids_checked

def generate_rsi_rso(number, inner, outer):
    """
    Generate arrays of sky annulus sizes

    Parameters
    ----------
    number : int
        The number of sky annuli required
    inner : int
        The inner radius of the sky annuli in pixels
    outer : int
        The outer radius of the sky annuli in pixels

    Returns
    -------
    rsi : array-like
        Array of inner sky annuli radii
    rso : array-like
        Array of outer sky annuli radii

    Raises
    ------
    None
    """
    rsi = np.ones(number)*inner
    rso = np.ones(number)*outer
    return rsi, rso

def source_extract(filename, sigma, rad_sky_inner, rad_sky_outer,
                   output=False, seg_map=False, check_image_boundary=False):
    """
    Measure the sky background and locate and extract all stars
    in the image supplied

    Parameters
    ----------
    filename : str
        The name of the file to run sep.extract on
    sigma : int, optional
        The number of sigma a detection must be above the background
        to be flagged as a star
    rad_sky_inner : int, optional
        The inner radius in pixels of the sky annulus
    rad_sky_outer : int, optional
        The outer radius in pixels of the sky annulus
    output : str, optional
        The name of the output file containing the positions of
        extracted stars
        Default = False (no output file)
    seg_map : bool, optional
        Toggle to select generating a SEP segmentation map of the
        image
        Default = False (no seg_map)
    check_image_boundary : bool
        Check the extracted positions are not too close to the edge
        Default = False (no boundary check)

    Returns
    -------
    object_ids : array-like
        List of object IDs
    x : array-like
        List of the X positions
    y : array-like
        List of the Y positions
    rsi : array-like
        List of inner sky annulus radii
    rso : array-like
        List of outer sky annulus radii
    segmentation_map : array-like, optional
        SEP segmentation map of image analysed

    Raises
    ------
    None
    """
    with fits.open(filename) as fitsfile:
        data = fitsfile[0].data.astype(np.float64)
        width = fitsfile[0].header['NAXIS1']
        height = fitsfile[0].header['NAXIS2']
        border = 50
    bkg = sep.Background(data)
    thresh = sigma * bkg.globalrms
    if not seg_map:
        objects = sep.extract(data-bkg, thresh)
        segmentation_map = None
    else:
        objects, segmentation_map = sep.extract(data-bkg, thresh, segmentation_map=seg_map)
    x = objects['x']
    y = objects['y']
    object_ids = list(np.arange(0, len(x)))

    # check the stars are not too close to the edge
    if check_image_boundary:
        x, y, object_ids = check_image_boundaries(x, y, object_ids, width, height, border)

    # generate sky annuli
    n_objects = len(x)
    rsi, rso = generate_rsi_rso(n_objects, rad_sky_inner, rad_sky_outer)

    if output:
        np.savetxt('{}.allstars'.format(output),
                   np.c_[x, y], fmt='%.2f  %.2f', header='x  y')
    return object_ids, x, y, rsi, rso, segmentation_map

def recenter_stars(xinit, yinit, source_x, source_y, sep_shift=1):
    """
    Take a set of XY measurements and compare them to source extracted
    positions to find better centroids. This is typically done if using
    a catalogue to set photometry apertures.

    Parameters
    ----------
    xinit : array-like
        List of initial X position guesses
    yinit : array-like
        List of initial Y position guesses
    source_x : array-like
        List of source extracted X positions to compare against
    source_y : array-like
        List of source extracted Y positions to compare against
    sep_shift : int, optional
        Maximium difference between guess and position found
        by SEP. If the difference > sep_shift, abort
        Default = 1

    Returns
    -------
    xo : array-like
        Updated list of X positions
    yo : array-like
        Updated list of Y positions

    Raises
    ------
    None
    """
    xo, yo = [], []
    for i, j in zip(xinit, yinit):
        diff_x = abs(source_x - i)
        diff_y = abs(source_y - j)
        radius = np.sqrt(diff_x**2+diff_y**2)
        match = np.where(radius == min(radius))[0][0]
        if diff_x[match] >= sep_shift or diff_y[match] >= sep_shift:
            print('Large shifts in recentroiding, check! skipping...')
            return [], []
        print(source_x[match], source_y[match])
        xo.append(source_x[match])
        yo.append(source_y[match])
    return np.array(xo), np.array(yo)

def get_light_travel_times(ra, dec, time_to_correct):
    """
    Get the light travel times to the helio- and
    barycentres

    Parameters
    ----------
    ra : str
        The Right Ascension of the target in hourangle
        e.g. 16:00:00
    dec : str
        The Declination of the target in degrees
        e.g. +20:00:00
    time_to_correct : astropy.Time object
        The time of observation to correct. The astropy.Time
        object must have been initialised with an EarthLocation

    Returns
    -------
    ltt_bary : float
        The light travel time to the barycentre
    ltt_helio : float
        The light travel time to the heliocentre

    Raises
    ------
    None
    """
    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = time_to_correct.light_travel_time(target)
    ltt_helio = time_to_correct.light_travel_time(target, 'heliocentric')
    return ltt_bary, ltt_helio
