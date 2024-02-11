"""
Housekeeping functions for reducing data
"""
import glob as g
from datetime import (
    datetime,
    timedelta
    )
from collections import defaultdict
import fitsio

def load_fits_image(filename, ext=0, force_float=False):
    """
    Read a fits image and header with fitsio

    Parameters
    ----------
    filename : str
        filename to load
    ext : int
        extension to load
    force_float : bool
        force image data to be float on load

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
    if force_float:
        data = data.astype(float)
    return data, header

def load_fits_table(filename, ext=0):
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
    if header:
        fitsio.write(filename, data, header=header, clobber=clobber)
    else:
        fitsio.write(filename, data, clobber=clobber)

def get_image_file_collection(instrument_config, directory, glob_exclude):
    """
    Make a dict of what images are present in a given
    directory. Replaces ccdproc.ImageFileCollection

    Parameters
    ----------
    instrument_config : dict
        Information about instrumentation configuration
    directory : str
        The path to the folder to analyse for images
    glob_exclude : str
        The wildcard str to exclude any files, typically master calibrations

    Returns
    -------
        image_file_collection : dict
        Dict of lists of filenames of varying types

    Raises
    ------
    None
    """
    # get a list of fits images in specific folder
    templist = sorted(g.glob(f"{directory}/*.f*"))
    # remove anything that was requested
    exclude_list = sorted(g.glob(f"{directory}/{glob_exclude}"))
    filtered_templist = [i for i in templist if i not in exclude_list]

    # storage for info
    biases = []
    darks = []
    flats = defaultdict(list)
    science = defaultdict(dict)

    # fetch all the header keyword info
    image_type_keyword = instrument_config['imager']['image_type_keyword']
    filter_keyword = instrument_config['imager']['filter_keyword']
    bias_keyword = instrument_config['imager']['bias_keyword']
    dark_keyword = instrument_config['imager']['dark_keyword']
    flat_keyword = instrument_config['imager']['flat_keyword']
    image_keyword = instrument_config['imager']['image_keyword']
    object_keyword = instrument_config['imager']['object_keyword']

    # loop over all images
    for image in filtered_templist:
        _, header = load_fits_image(image)
        image_type = header[image_type_keyword]

        # keep a list of biases
        if image_type == bias_keyword:
            biases.append(image)
        # keep a list of darks
        elif image_type == dark_keyword:
            darks.append(image)
        # keep a nested list of all flats
        elif image_type == flat_keyword:
            filt = header[filter_keyword]
            flats[filt].append(image)
        # keep a nested list of science/filter objects
        elif image_type == image_keyword:
            obj = header[object_keyword]
            if obj not in science:
                science[obj] = defaultdict(list)
            filt = header[filter_keyword]
            science[obj][filt].append(image)

    # build the final image file collection dict
    image_file_collection = {'biases': biases,
                             'darks': darks,
                             'flats': flats,
                             'science': science}
    return image_file_collection

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
    _, hdr = load_fits_image(filename)
    return hdr[object_keyword]

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
    _, hdr = load_fits_image(filename)
    date_obs = datetime.strptime(hdr[dateobs_keyword].split('.')[0],
                                 "%Y-%m-%dT%H:%M:%S")
    if date_obs.hour < 12:
        date_obs = date_obs - timedelta(days=1)
    return date_obs.strftime('%Y%m%d')
