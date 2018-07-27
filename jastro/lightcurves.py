"""
A series of functions to handle light curves
"""
import numpy as np

def lc_mags_to_flux(mags, mags_err):
    """
    Take 2 arrays, light curve and errors
    and convert them from differential magnitudes
    back to relative fluxes

    Rolling back these equations:
        mags = - 2.5 * log10(flux)
        mag_err = (2.5/log(10))*(flux_err/flux)
    """
    flux = 10.**(mags / -2.5)
    flux_err = (mags_err/(2.5/np.log(10))) * flux
    return flux, flux_err
