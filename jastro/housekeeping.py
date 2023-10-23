"""
Housekeeping functions for reducing data
"""
import ccdproc

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
