"""
Housekeeping functions for reducing data
"""
from datetime import (
    datetime,
    timedelta
    )
from astropy.io import fits
import fitsio
import ccdproc


def load_fits_image(filename, ext=0):
    """
    Read a fits image and header with fitsio

    Parameters
    ----------
    filename : str
        filename to load
    ext : int
        extension to load

    Returns
    -------
    data : array
        data from the corresponding extension
    header : fitsio.header.FITSHDR
        list from file header

    Raises
    ------
    None
    """
    data, header = fitsio.read(filename, header=True, ext=ext)
    return data, header

def write_fits_image(filename, data, header, clobber=True):
    """
    Write a fits image to disc with fitsio

    Parameters
    ----------
    filename : str
        name | path of image to write
    data : array
        data to save
    header : list of dicts
        name, value, comment for each entry
    clobber : boolean
        overwrite existing image?
        default = True

    Returns
    -------
    None

    Raises
    ------
    None
    """
    fitsio.write(filename, data, header=header, clobber=clobber)

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

def get_target_id(filename, object_keyword='OBJECT'):
    """
    Get the target ID from the reference image

    Parameters
    ----------
    filename : str
        Name of the image from which to determine the night ID from
    object_keyword : str, optional
        Header keyword for the target ID
        Default = 'OBJECT'

    Returns
    -------
    target_id : str
        Name of the target from the reference image

    Raises
    ------
    None
    """
    with fits.open(filename) as fitsfile:
        target_id = fitsfile[0].header[object_keyword]
    return target_id

def get_night_id(filename, dateobs_keyword='DATE-OBS'):
    """
    Get the YYYYMMDD night ID from the reference image

    Parameters
    ----------
    filename : str
        Name of the image from which to determine the night ID from
    dateobs_keyword : str, optional
        Header keyword for the date of observation
        Default = 'DATE-OBS'

    Returns
    -------
    night : str
        YYYYMMDD night of observation

    Raises
    ------
    None
    """
    with fits.open(filename) as fitsfile:
        date_obs = datetime.strptime(fitsfile[0].header[dateobs_keyword].split('.')[0], "%Y-%m-%dT%H:%M:%S")
    if date_obs.hour < 12:
        date_obs = date_obs - timedelta(days=1)
    return date_obs.strftime('%Y%m%d')
