"""
Functions for loading config info
"""
import json

# pylint: disable=invalid-name
# pylint: disable=unspecified-encoding

def load(filename):
    """
    Load a json file

    Parameters
    ----------
    filename : str
        name of file to load

    Returns
    -------
    config : dict
        dict of json file contents

    Raises
    ------
    None
    """
    with open(filename, "r") as ff:
        config = json.load(ff)
    return config
