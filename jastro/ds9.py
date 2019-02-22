"""
DS9 functionality
"""
import os
import re
import time
import platform
from collections import defaultdict
import numpy as np
import pyregion

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# pylint: disable = invalid-name
# pylint: disable = line-too-long
# pylint: disable = too-many-arguments
# pylint: disable = too-many-locals

def check_for_ds9(ds9_name):
    """
    Find DS9 window with ds9_name if it exists

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to check for

    Returns
    -------
    ds9 : bool
        Boolean showing if DS9 was found

    Raises
    ------
    None
    """
    ds9 = False
    if platform.system().lower() == 'darwin':
        processes = os.popen('ps aux | grep {}'.format(ds9_name)).readlines()
    else:
        processes = os.popen('ps aux --cols 1024 | grep {}'.format(ds9_name)).readlines()
    for process in processes:
        if "ds9" in process and " -title {}".format(ds9_name) in process:
            print('Hurray! Found DS9 window process')
            print('Waiting 20s to be sure ds9 is open...')
            time.sleep(20)
            ds9 = True
            break
    else:
        print('Boo! Where the hell is the DS9 window?')
    return ds9

def start_ds9(ds9_name):
    """
    Start up a DS9 window. Doing it like this is much
    more robust than using pyds9.

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to open

    Returns
    -------
    None

    Raises
    ------
    None
    """
    os.system('ds9 -title {} &'.format(ds9_name))

def display_ds9(ds9_name, image_id):
    """
    Display an image in DS9

    Parameters
    ----------
    ds9_name : str
        ID of the ds9 window to display an image in
    image_id : str
        Name of the FITS image to display

    Returns
    -------
    None

    Raises
    ------
    None
    """
    os.system('xpaset {} fits < {}'.format(ds9_name, image_id))

def set_ds9(ds9_name, command):
    """
    Set command for updating a ds9 window

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to update
    command : str
        xpaset command for updating DS9 window
        see http://ds9.si.edu/doc/ref/xpa.html for more details
    Returns
    -------
    None

    Raises
    ------
    None
    """
    os.system('xpaset -p {} {}'.format(ds9_name, command))

def setup_ds9(ds9_name, filename=None):
    """
    Configure the DS9 image viewing window

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to set up
    filename : str
        Name of the file to display in DS9

    Returns
    -------
    None

    Raises
    ------
    None
    """
    print('Checking for DS9...')
    if not check_for_ds9(ds9_name):
        print('No DS9 found, starting DS9...')
        start_ds9(ds9_name)
        while not check_for_ds9(ds9_name):
            print('Checking for DS9...')
            time.sleep(1)

    time.sleep(5)
    set_ds9(ds9_name, 'scale mode zscale')
    set_ds9(ds9_name, 'preserve scale yes')
    set_ds9(ds9_name, 'preserve pan yes')
    if filename:
        display_ds9(ds9_name, filename)
        set_ds9(ds9_name, 'zoom to fit')
        region_file = "{0:s}.reg".format(filename.split('.')[0])
        if os.path.exists(region_file):
            set_ds9(ds9_name, 'regions < {0:s}'.format(region_file))

def output_ds9_region_file(reference_image_prefix, object_ids, x, y,
                           r_aper, rsi=None, rso=None, index_offset=1.0):
    """
    Take the output of convertWorldToPix and make a region file
    This region file can then be overlaid and edited by hand
    (mainly the sky annuli if not provided)
    Parameters
    ----------
    reference_image_prefix : str
        Prefix of the region file name (the part you want to
        come before .reg)
    object_ids : array-like
        List of object IDs
    x : array-like
        X positions of the stars
    y : array-like
        Y positions of the stars
    r_aper : int
        Radius of the photometry aperture
    rsi : array-like, optional
        Radii of inner edge of sky annuli
        Default = None
    rso : array-like, optional
        Radii of outer edge of sky annuli
        Default = None

    Returns
    -------
    None

    Raises
    ------
    None
    """
    outfile = '{}.reg'.format(reference_image_prefix)
    if os.path.exists(outfile):
        print('WARNING: CLOBBERING {}'.format(outfile))
    f = open(outfile, 'w')
    region_header = ('# Region file format: DS9 version 4.1\n'
                     'global color=red dashlist=8 3 width=1 '
                     'font="helvetica 10 normal roman" select=1 '
                     'highlite=1 dash=0 fixed=0 edit=1 move=1 '
                     'delete=1 include=1 source=1\nphysical\n')
    f.write(region_header)
    if not rso or not rsi:
        rsi = np.array([20]*len(x))
        rso = np.array([25]*len(x))
    for i, j, k, l, m in zip(x, y, rsi, rso, object_ids):
        line = 'annulus({},{},{},{}) # text={{{}}}\n'.format(i+index_offset,
                                                             j+index_offset,
                                                             k, l, m)
        line2 = 'circle({},{},{}) \n'.format(i+index_offset,
                                             j+index_offset,
                                             r_aper)
        f.write(line)
        f.write(line2)
    f.close()

def read_ds9_region_file(filename, fits_file_extension='.f', index_offset=1.0):
    """
    Read in ds9 regions from the ref image region file
    (if it exists)

    We assume that targets have proper named regions
    and that comparison stars have integer IDs
    (e.g. Comps = 1, 2, 3, 4 and Target = WASP-148b)

    Parameters
    ----------
    filename : str
        The name of the file to read regions for
    fits_file_extension : str, optional
        The filename extension
        e.g. '.fit', '.fits' or '.fts'
        Default = '.fts'

    Returns
    -------
    x : array-like
        List of X positions
    y : array-like
        List of Y positions
    rsi : array-like
        List of inner radii of sky annuli
    rso : array-like
        List of outer radii of sky annuli

    Raises
    ------
    None
    """
    filename = '{0:s}.reg'.format(filename.split(fits_file_extension)[0])
    print('Looking for {0:s}'.format(filename))
    try:
        f = pyregion.open(filename)
    except FileNotFoundError:
        print('File not found, do photometry of all stars')
        return [], [], [], []

    comparisons = defaultdict(list)
    x, y, rsi, rso = [], [], [], []
    region_regex = r'text={(?P<label>\w+)}'
    r = re.compile(region_regex)
    for i in f:
        match = r.match(i.comment)
        if match:
            print('Found star {0:s} in {1:s}'.format(match.group('label'), filename))
            try:
                comparisons[int(match.group('label'))] = ['{0:.2f}'.format(float(i.coord_list[0])-index_offset),
                                                          '{0:.2f}'.format(float(i.coord_list[1])-index_offset),
                                                          '{0:.2f}'.format(float(i.coord_list[2])),
                                                          '{0:.2f}'.format(float(i.coord_list[3]))]
            except ValueError:
                target = ['{0:.2f}'.format(float(i.coord_list[0])-index_offset),
                          '{0:.2f}'.format(float(i.coord_list[1])-index_offset),
                          '{0:.2f}'.format(float(i.coord_list[2])),
                          '{0:.2f}'.format(float(i.coord_list[3]))]
    # comparisons
    for i in range(0, len(comparisons)):
        x.append(float(comparisons[i+1][0]))
        y.append(float(comparisons[i+1][1]))
        rsi.append(float(comparisons[i+1][2]))
        rso.append(float(comparisons[i+1][3]))
    # target
    x.append(float(target[0]))
    y.append(float(target[1]))
    rsi.append(float(target[2]))
    rso.append(float(target[3]))
    return np.array(x), np.array(y), np.array(rsi), np.array(rso)

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
    return round(data[y-radius:y+radius, x-radius:x+radius].ravel().max(), 2)

