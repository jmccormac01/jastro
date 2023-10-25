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

# pylint: disable = invalid-name
# pylint: disable = line-too-long
# pylint: disable = too-many-arguments
# pylint: disable = too-many-locals

def is_running(ds9_name):
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
        processes = os.popen(f'ps aux | grep {ds9_name}').readlines()
    else:
        processes = os.popen(f'ps aux --cols 1024 | grep {ds9_name}').readlines()
    for process in processes:
        if "ds9" in process and f" -title {ds9_name}" in process:
            print('Hurray! Found DS9 window process')
            print('Waiting 5s to be sure ds9 is open...')
            time.sleep(5)
            ds9 = True
            break
    else:
        print('Boo! Where the hell is the DS9 window?')
    return ds9

def start(ds9_name):
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
    os.system(f'/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title {ds9_name} &')

def display(ds9_name, image_id):
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
    os.system(f'xpaset {ds9_name} fits < {image_id}')

def dset(ds9_name, command):
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
    os.system(f'xpaset -p {ds9_name} {command}')

def setup(ds9_name, filename=None):
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
    if not is_running(ds9_name):
        print('No DS9 found, starting DS9...')
        start(ds9_name)
        while not is_running(ds9_name):
            print('Checking for DS9...')
            time.sleep(1)

    time.sleep(5)
    dset(ds9_name, 'scale mode zscale')
    dset(ds9_name, 'preserve scale yes')
    dset(ds9_name, 'preserve pan yes')
    if filename:
        display(ds9_name, filename)
        dset(ds9_name, 'zoom to fit')
        prefix = filename.split('.')[0]
        region_file = f"{prefix}.reg"
        if os.path.exists(region_file):
            dset(ds9_name, f'regions < {region_file}')

def output_region_file(reference_image_prefix, object_ids, x, y,
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
    outfile = f'{reference_image_prefix}.reg'
    if os.path.exists(outfile):
        print(f'WARNING: CLOBBERING {outfile}')
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
        line = f'annulus({i+index_offset},{j+index_offset},{k},{l}) # text={{{m}}}\n'
        line2 = f'circle({i+index_offset},{j+index_offset},{r_aper}) \n'
        f.write(line)
        f.write(line2)
    f.close()

def read_region_file(region_filename, index_offset=1):
    """
    Read in ds9 regions from the ref image region file
    (if it exists)

    We assume that targets have proper named regions
    and that comparison stars have integer IDs
    (e.g. Comps = 1, 2, 3, 4 and Target = WASP-148b)

    Parameters
    ----------
    region_filename : str
        file name with sky annulii
    index_offsex : int
        difference between python and ds9 indexing

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
    print(f'Looking for {region_filename}')
    try:
        f = pyregion.open(region_filename)
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
            lab = match.group('label')
            print(f'Found star {lab} in {region_filename}')
            try:
                comparisons[int(lab)] = ['{0:.2f}'.format(float(i.coord_list[0])-index_offset),
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
