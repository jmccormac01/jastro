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

def phase_times(times, epoch, period, phase_offset=0.0):
    """
    Take a list of times, an epoch, a period and
    convert the times into phase space.

    An optional offset can be supplied to shift phase 0
    """
    return (((times - epoch)/period)+phase_offset)%1

def pc_bin(time, mmi, bin_width, clip_empty_bins=True):
    bin_edges = np.arange(np.min(time), np.max(time), bin_width)
    digitized = np.digitize(time, bin_edges)
    binned_time = (bin_edges[1:] + bin_edges[:-1]) / 2
    binned_mmi = np.array([np.nan if len(mmi[digitized == i]) == 0 else mmi[digitized == i].mean() for i in range(1, len(bin_edges))])
    if clip_empty_bins:
        binned_time = binned_time[~np.isnan(binned_mmi)]
        binned_mmi = binned_mmi[~np.isnan(binned_mmi)]
    return (binned_time, binned_mmi)
