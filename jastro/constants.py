"""
List of constants used frequently in astronomical calculations
"""

class Constants(object):
    """
    Common suite of useful constants in astronomy

    Solar mass, radius and luminosity values are from
    Harmanec & Prsa arXiv:1106.1508v2
    """
    G = 6.67408E-11 # m3 kg^-1 s^-2
    Msun = 1.988547E30 # kg
    Rsun = 6.95508E8 # m
    Lsun = 3.846E26 # W
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
