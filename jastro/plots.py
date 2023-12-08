"""
Functions for plotting aperture photometry
"""
import matplotlib.pyplot as plt
import numpy as np

# pylint: disable = invalid-name
# pylint: disable = redefined-outer-name
# pylint: disable = no-member
# pylint: disable = too-many-locals
# pylint: disable = too-many-arguments
# pylint: disable = unused-variable
# pylint: disable = line-too-long

def plot_max_pixel_values(t, comparisons_max_pix, target_max_pix,
                          target_id, filt, aperture_radius):
    """
    Plot the maximum value in the box around the photometry
    aperture for each comparison and the target.
    This can be used to flag bad comparison stars
    """
    fig, ax = plt.subplots(1, figsize=(10, 10))
    leg = []
    for i in range(0, len(comparisons_max_pix)):
        ax.plot(t, comparisons_max_pix[i], '.')
        leg.append(i+1)
    leg.append('Target')
    ax.plot(t, target_max_pix, 'o')
    ax.set_xlabel('JD')
    ax.set_ylabel('Reduced Flux')
    ax.set_title('Maximum Pixel Values')
    ax.legend((leg), loc=1)
    fig.savefig(f'{target_id.replace(" ", "-")}_{filt}_A{aperture_radius}_max_pixel_values.png')
    plt.show()

def plot_comparison_stars(t, comparisons, target_id, filt, aperture_radius):
    """
    Compare the comparison stars to each other
    to look for variable stars. Fit a 2nd order
    polynomial to them to crudely flatten airmass
    and allow measurement of the RMS
    """
    # set up the plotting grid geometry
    NCOLS = 2
    COL_WIDTH = 2
    PLOT_COLS = NCOLS * COL_WIDTH
    SUB_ROWS = 1
    ROW_WIDTH = 1
    NROWS = int(np.ceil(len(comparisons)/float(NCOLS))*SUB_ROWS)
    PLOT_ROWS = NROWS * ROW_WIDTH

    # arrays to hold data
    coeffs = np.empty((len(comparisons), 3))
    cbesty = np.empty((len(comparisons), len(comparisons[0])))
    crms = np.empty(len(comparisons))

    for i in range(0, len(comparisons)):
        coeff = np.polyfit(t, comparisons[i], 2)
        coeffs[i] = np.array(coeff)
        cbesty[i] = np.polyval(coeffs[i], t)
        # flatten and get RMS
        crms[i] = np.std(comparisons[i] / cbesty[i])
    crms_avg = np.average(crms)
    print(f'Average RMS: {crms_avg:.4f}')

    # now do the plotting
    for i in range(0, len(comparisons)):
        fig = plt.figure(i+1, figsize=(10, 10))
        fig.suptitle(f'Comparison Star {i+1:d}')
        c = 0
        r = 0
        for j in range(0, len(comparisons)):
            ax = plt.subplot2grid((PLOT_ROWS, PLOT_COLS), (r, c*COL_WIDTH),
                                  colspan=COL_WIDTH, rowspan=ROW_WIDTH)
            ratio = comparisons[i]/comparisons[j]
            ratio_normalised = ratio / np.average(ratio[:10])
            ax.plot(t, ratio_normalised, 'r.')
            ax.set_title(f'ratio with comparison {j+1:d}')
            c += 1
            if c > NCOLS - 1:
                c = 0
                r += SUB_ROWS
        fig.savefig(f'{target_id.replace(" ", "-")}_{filt}_A{aperture_radius}_comparisons_wrt_{i+1}.png')
    plt.show()

def plot_star_fluxes(t, comparisons, target, target_id, filt, aperture_radius):
    """
    Plot the raw fluxes for each object
    """
    fig, ax = plt.subplots(1, figsize=(10, 10))
    leg = []
    for i in range(0, len(comparisons)):
        ax.plot(t, comparisons[i], '.')
        leg.append(i+1)
    leg.append('Target')
    ax.plot(t, target, 'o')
    ax.set_xlabel('JD')
    ax.set_ylabel('Flux')
    ax.set_title('Raw Fluxes')
    ax.legend((leg), loc=1)
    fig.savefig(f'{target_id.replace(" ", "-")}_{filt}_A{aperture_radius}_raw_fluxes.png')
    plt.show()
