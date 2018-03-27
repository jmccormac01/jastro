"""
List of constants used frequently in astronomical calculations
"""

class Constants(object):
    """
    Common suite of useful constants in astronomy
    """
    G = 6.67408E-11 # m3 kg^-1 s^-2
    Msun = 1.989E30 # kg
    Rsun = 695700000 # m
    Mjup = 1.898E27 # kg
    Rjup = 69911000 # m
    Mproton = 1.6726219E-27 # kg
    boltzmann = 1.38064852E-23 # m^2 kg s^-2 K^-1
    AU = 1.496E11 # m
    sec_d = 24*60*60 # s

    def _list_constants():
        """
        Print a list of all the constants above
        """
        c = Constants
        print(" ".join([i for i in vars(c) if not i.startswith("_")]))
