"""
A series of functions to handle light curves
"""
from collections import (
    OrderedDict,
    defaultdict
    )
import numpy as np

# pylint: disable=invalid-name

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

def lc_flux_to_mags(flux, flux_err):
    """
    Take 2 arrays, light curve and errors
    and convert them from differential magnitudes
    back to relative fluxes

    Applying these equations:
        mags = - 2.5 * log10(flux)
        mag_err = (2.5/log(10))*(flux_err/flux)
    """
    mags = -2.5*np.log10(flux)
    mags_err = (2.5/np.log(10))*(flux_err/flux)
    return mags, mags_err

def phase_times(times, epoch, period, phase_offset=0.0):
    """
    Take a list of times, an epoch, a period and
    convert the times into phase space.

    An optional offset can be supplied to shift phase 0
    """
    return (((times - epoch)/period)+phase_offset)%1

def pc_bin(time, flux, error, bin_width, clip_empty_bins=True, mode="mean"):
    """
    """
    bin_edges = np.arange(np.min(time), np.max(time), bin_width)
    digitized = np.digitize(time, bin_edges)
    binned_time = (bin_edges[1:] + bin_edges[:-1]) / 2

    if mode == "median":
        binned_flux = np.array([np.nan if len(flux[digitized == i]) == 0 else np.median(flux[digitized == i]) for i in range(1, len(bin_edges))])
        binned_error = np.array([np.nan if len(error[digitized == i]) == 0 else np.sqrt(np.sum(error[digitized == i]**2.))/len(error[digitized == i]) for i in range(1, len(bin_edges))])
    else:
        binned_flux = np.array([np.nan if len(flux[digitized == i]) == 0 else flux[digitized == i].mean() for i in range(1, len(bin_edges))])
        binned_error = np.array([np.nan if len(error[digitized == i]) == 0 else np.sqrt(np.sum(error[digitized == i]**2.))/len(error[digitized == i]) for i in range(1, len(bin_edges))])

    if clip_empty_bins:
        binned_time = binned_time[~np.isnan(binned_flux)]
        binned_flux = binned_flux[~np.isnan(binned_flux)]
        binned_error = binned_error[~np.isnan(binned_flux)]
    return (binned_time, binned_flux, binned_error)

def extract_nights_with_transits(times, flux, err, epoch,
                                 period, t14, transit_type='full'):
    """
    Take in a time series and ephmeris and
    return only the nights with transits
    """
    # split the data into separate transits
    nights = OrderedDict()
    skip, hit = 0, 0
    phase_lim = (t14/2.)/period
    for i, t in enumerate(times):
        ph = (t - epoch)/period
        if ph%1 >= 0.90 or ph%1 <= 0.1:
            try:
                nights[round(ph)]['hjd'].append(t)
                nights[round(ph)]['flux'].append(flux[i])
                nights[round(ph)]['phase'].append(ph%1)
                nights[round(ph)]['error'].append(err[i])
            except KeyError:
                nights[round(ph)] = defaultdict(list)
                nights[round(ph)]['hjd'].append(t)
                nights[round(ph)]['flux'].append(flux[i])
                nights[round(ph)]['phase'].append(ph%1)
                nights[round(ph)]['error'].append(err[i])
            hit += 1
        else:
            skip += 1

    if transit_type == 'partial':
        # Do a second pass to see if none of the data points are in transit
        print('Checking full and partial transits...')
        to_delete = []
        for night in nights:
            hits = 0
            for p in nights[night]['phase']:
                if p >= 1-phase_lim or p <= phase_lim:
                    hits += 1
            if hits == 0:
                print("Delete transit {}, no points within phase_lim".format(night))
                to_delete.append(night)
    elif transit_type == 'full':
        print('Checking for only full transits...')
        to_delete = []
        for night in nights:
            hits_l = 0
            hits_p = 0
            for p in nights[night]['phase']:
                if p <= 1-phase_lim and p > 0.5:
                    hits_l += 1
                elif p >= phase_lim and p < 0.5:
                    hits_p += 1
            if hits_l < 10 or hits_p < 10:
                to_delete.append(night)
    else:
        return None

    # do the final cut
    print("N nights pre-clean-up: {}".format(len(nights)))
    for d in to_delete:
        del nights[d]
    print("N nights post-clean-up: {}".format(len(nights)))
    return nights
