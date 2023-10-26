"""
A series of functions to handle light curves
"""
from collections import (
    OrderedDict,
    defaultdict
    )
import numpy as np
import matplotlib.pyplot as plt

# pylint: disable=invalid-name

def mags_to_flux(mags, mags_err):
    """
    Convert diff mags to relative fluxes

    Rolling back these equations:
        mags = - 2.5 * log10(flux)
        mag_err = (2.5/log(10))*(flux_err/flux)

    Parameters
    ----------
    mags : array
        diff mags to convert
    mags_err : array
        error on diff mags to convert

    Returns
    -------
    flux : array
        relative fluxes from diff mags
    flux_err : array
        errors on relative fluxes from diff mags

    Raises
    ------
    None
    """
    flux = 10.**(mags / -2.5)
    flux_err = (mags_err/(2.5/np.log(10))) * flux
    return flux, flux_err

def flux_to_mags(flux, flux_err):
    """
    Convert relative fluxes to diff mags

    Applying these equations:
        mags = - 2.5 * log10(flux)
        mag_err = (2.5/log(10))*(flux_err/flux)

    Parameters
    ----------
    flux : array
        relative fluxes from diff mags
    flux_err : array
        errors on relative fluxes from diff mags

    Returns
    -------
    mags : array
        diff mags to convert
    mags_err : array
        error on diff mags to convert

    Raises
    ------
    None
    """
    mags = -2.5*np.log10(flux)
    mags_err = (2.5/np.log(10))*(flux_err/flux)
    return mags, mags_err

def phase_times(times, epoch, period, phase_offset=0.0):
    """
    Take a list of times, an epoch, a period and
    convert the times into phase space.

    An optional offset can be supplied to shift phase 0

    Parameters
    ----------
    times : array
        Array of times to phase
    epoch : floar
        T0 for phasing
    period : float
        Period to fold on
    phase_offset : float
        Shift phase by fraction of a period

    Returns
    -------
    phase : array
        phased time array

    Raises
    ------
    None
    """
    return (((times - epoch)/period)+phase_offset)%1

def bin_time_flux_error(time, flux, error, bin_fact):
    """
    Use reshape to bin light curve data, clip under filled bins
    Works with 2D arrays of flux and errors

    Note: under filled bins are clipped off the end of the series

    Parameters
    ----------
    times : array
        Array of times to bin
    flux : array
        Array of flux values to bin
    error : array
        Array of error values to bin
    bin_fact : int
        Number of measurements to combine

    Returns
    -------
    times_b : array
        Binned times
    flux_b : array
        Binned fluxes
    error_b : array
        Binned errors

    Raises
    ------
    None
    """
    n_binned = int(len(time)/bin_fact)
    clip = n_binned*bin_fact
    time_b = np.average(time[:clip].reshape(n_binned, bin_fact), axis=1)
    # determine if 1 or 2d flux/err inputs
    if len(flux.shape) == 1:
        flux_b = np.average(flux[:clip].reshape(n_binned, bin_fact), axis=1)
        error_b = np.average(error[:clip].reshape(n_binned, bin_fact), axis=1)
    else:
        # assumed 2d with 1 row per star
        n_stars = len(flux)
        flux_b = np.average(flux[:clip].reshape((n_stars, n_binned, bin_fact)), axis=2)
        error_b = np.average(error[:clip].reshape((n_stars, n_binned, bin_fact)), axis=2)
    return time_b, flux_b, error_b

def extract_nights_with_transits(times, flux, err, epoch,
                                 period, t14, transit_type='full'):
    """
    Take in a time series and ephemeris and
    return only the nights with transits

    Parameters
    ----------
    times : array
        Array of times
    flux : array
        Array of flux values
    err : array
        Array of error values
    epoch : floar
        T0 for extracting ephem
    period : float
        Period for extracting ephem
    t14 : float
        Transit duration for extracting ephem
    transit_type : str
        Extract 'full' or 'partial' transits?

    Returns
    -------
    nights : dict
        Collection of lcs from nights with transits

    Raises
    ------
    None
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
                print(f"Delete transit {night}, no points within phase_lim")
                to_delete.append(night)
    elif transit_type == 'full':
        print('Checking for only full transits...')
        to_delete = []
        for night in nights:
            hits_l = 0
            hits_p = 0
            for p in nights[night]['phase']:
                if 0.5 < p <= 1-phase_lim:
                    hits_l += 1
                elif phase_lim <= p < 0.5:
                    hits_p += 1
            if hits_l < 10 or hits_p < 10:
                to_delete.append(night)
    else:
        return None

    # do the final cut
    print(f"N nights pre-clean-up: {len(nights)}")
    for d in to_delete:
        del nights[d]
    print(f"N nights post-clean-up: {len(nights)}")
    return nights

def normalise(filt, t, t0, lightcurve, lightcurve_err, r_aper, bin_fact,
              target_id, night_id, fit_type=1, fit_low=None,
              fit_high=None, ylim=None):
    """
    Basic fitting of transit lightcurve until we have
    something more sophisticated. Not happy with this
    but it will do for now (If I am reading this 6 months
    from now - 20160708, tut tut!)

    Parameters
    ----------
    filt : str
        Filter used for observations
    t : array-like
        Times of data points
    t0 : int
        Integer day of the first data point
    lightcurve : array-like
        Lightcurve data array
    lightcurve_err : array-like
        Error on Lightcurve
    r_aper : float
        Radius of photometry aperture
    bin_fact : int
        Number of data points to bin together
    target_id : str
        Name of object
    night_id : str
        Night of observation
    fit_type : int, optional
        Order of the polynomal fit to the out of transit data
        Default = 1 (Only 1 or 2 supported)
    fit_low : float
        Time below which to fit during normalisation
        Time assumed in fractional day
        Default = None
    fit_high : float
        Time above which to fit during normalisation
        Time assumed in fractional day
        Default = None
    ylim : str | None, optional
        Comma-separated string giving the lower and upper
        bounds of the Y-axis in the final transit plot
        Default = None

    Returns
    -------
    lightcurve_n : array-like
        Normalised lightcurve, still in flux units
    lightcurve_err_n : array-like
        Error on lightcurve_n

    Raises
    ------
    None
    """
    # determine which sections of the data to fit
    if fit_low is not None and fit_high is not None:
        print('Fitting pre and post...')
        index = np.where(((t < fit_low) | (t > fit_high)))[0]
    elif fit_low is not None and fit_high is None:
        print('Fitting pre...')
        index = np.where(t < fit_low)[0]
    elif fit_low is None and fit_high is not None:
        print('Fitting post...')
        index = np.where(t > fit_high)[0]
    else:
        print('Fitting everything...')
        index = np.where(t > 0.0)[0]

    # crude fit
    if fit_type > 0:
        # REPLACE THIS WITH SOMETHING MUCH BETTER!
        # 1D or 2D fit to OOT data
        coeffs = np.polyfit(t[index], lightcurve[index], fit_type)
        besty = np.polyval(coeffs, t)
    elif fit_type == 0:
        # use an average of the fit_sect area to normalise
        # this is used for cases where only one OOT exists
        # and fitting a line cocks things up
        besty = np.full(len(t), np.average(lightcurve[index]))
    else:
        # use this to normalise fluxes from a ratio of previous night
        scale = input('Enter scaling factor: ')
        besty = np.full(len(t), float(scale))
    lightcurve_n = lightcurve / besty
    lightcurve_err_n = lightcurve_err / lightcurve_n
    rms = np.std(lightcurve_n[index])
    print(f'RMS-{fit_type}: {rms:.4f}')

    # bin up the data
    tb, lightcurve_nb, _ = bin_time_flux_error(t, lightcurve_n,
                                               lightcurve_err_n, bin_fact)
    # work out the binned indexes
    if fit_low is not None and fit_high is not None:
        index_b = np.where(((tb < fit_low) | (tb > fit_high)))[0]
    elif fit_low is not None and fit_high is None:
        index_b = np.where(tb < fit_low)[0]
    elif fit_low is None and fit_high is not None:
        index_b = np.where(tb > fit_high)[0]
    else:
        index_b = np.where(tb > 0.0)[0]

    # get the binned RMS
    rmsb = np.std(lightcurve_nb[index_b])
    print(f'RMS-{fit_type}b: {rmsb:.4f}')

    # define the ylim if required
    if ylim:
        y_llim, y_ulim = ylim.split(',')
        y_llim = float(y_llim)
        y_ulim = float(y_ulim)
    fig, ax = plt.subplots(3, sharex=True, figsize=(10, 10))
    # raw data
    ax[0].errorbar(t, lightcurve, yerr=lightcurve_err, fmt='.', color='r', ecolor='lightgrey')
    ax[0].plot(t, besty, 'k--')
    ax[0].set_ylabel('Target / Comparison')
    ax[0].legend((f'{fit_type}D polyfit', 'Data',), loc='best')
    ax[0].set_title('Raw Lightcurve')
    ax[0].set_xlim(min(t)-0.05, max(t)+0.07)
    # unbinned data
    ax[1].errorbar(t, lightcurve_n, yerr=lightcurve_err_n, fmt='.', color='r', ecolor='lightgrey')
    ax[1].set_ylabel('Normalised Flux')
    ax[1].set_title(f'Raw Lightcurve / {fit_type}D Model')
    ax[1].legend((f'RMS-{fit_type} = {rms:.4f}',), loc='best')
    if ylim:
        ax[1].set_ylim(y_llim, y_ulim)
    ax[1].set_xlim(min(t)-0.05, max(t)+0.07)
    # binned data
    ax[2].plot(tb, lightcurve_nb, 'r.')
    ax[2].set_ylabel('Normalised Flux (binned)')
    ax[2].set_xlabel(f'JD - {t0:d}+')
    ax[2].set_title(f'Raw Lightcurve / {fit_type}D Model (binned x {bin_fact:d})')
    ax[2].legend((f'RMS-{fit_type}b = {rmsb:.4f}',), loc='best')
    if ylim:
        ax[2].set_ylim(y_llim, y_ulim)
    ax[2].set_xlim(min(t)-0.05, max(t)+0.07)

    plotname = f"{target_id}-{filt}-{night_id}-F{fit_type}-A{r_aper}.png"
    fig.tight_layout()
    fig.savefig(plotname)
    plt.show()
    return lightcurve_n, lightcurve_err_n
