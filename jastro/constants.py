"""
List of constants used frequently in astronomical calculations
"""
import math

# pylint: disable=invalid-name

class Constants(object):
    """
    Common suite of useful constants in astronomy

    Radii, masses etc for solar system bodies are taken
    from Mamajek et al. 2015
    https://arxiv.org/pdf/1510.07674.pdf
    """
    # the 2014 CODATA value is G = 6.67408 (±0.00031) × 10−11 m3 kg−1 s−2
    G = 6.67408E-11 # m3 kg^-1 s^-2
    # Sun
    GMsun = 1.3271244E20 # m^3 s^-2
    Msun = GMsun/G # kg
    Rsun = 6.957E8 # m
    Lsun = 3.828E26 # W
    Ssun = 1361. # W m^-2
    Teffsun = 5772. # K
    Rhosun = Msun / ((4./3.)*math.pi*(Rsun**3.))
    # Jupiter
    GMjup = 1.2668653E17 # m^3 s^-2
    Mjup = GMjup/G # kg
    Rpjup = 6.6854E7 # m - polar radius
    Rejup = 7.1492E7 # m - equatorial radius
    Rjup = Rejup # if not specified, we use the equatorial
    # Earth
    GMe = 3.986004E14 # m^3 s^-2
    Me = GMe/G # kg
    Rpe = 6.3568E6 # m - polar radius
    Ree = 6.3781E6 # m - equatorial radius
    Re = Ree # if not specified, we use the equatorial
    # Misc astronomy
    AU = 1.496E11 # m
    parsec = 3.0857E16 # m
    # Time
    sec_d = 24*60*60 # s
    # Physics
    Mproton = 1.6726219E-27 # kg
    stefboltz = 5.67E-8 # W m^-2 K^-4

    def _list_constants():
        """
        Print a list of all the constants above
        """
        c = Constants
        print(" ".join([i for i in vars(c) if not i.startswith("_")]))
