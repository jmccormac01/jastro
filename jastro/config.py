"""
Functions for loading config info
"""
import json

def load(filename):
    with open(filename, "r") as ff:
        config = json.load(ff)
    return config
