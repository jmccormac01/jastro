"""
Functions for extracting photometry
"""
import numpy as np
import sep
from jastro.ds9 import dset
from jastro.reduce import find_max_pixel_value

# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=c-extension-no-member
# pylint: disable=line-too-long
# pylint: disable=consider-using-f-string

# TODO: redo the output format later

def phot(data, shift, x, y, rsi, rso, aper_radii, filename, jd, bjd, hjd,
         phot_filename_prefix="rtp", ds9_name=None, draw_regions=False,
         gain=1.00, index_offset=1):
    """
    Measure the photometry of the given image at the positions
    requested. Output the results to a photometry file and
    display the apertures in DS9 if required.

    Parameters
    ----------
    data : ccdproc.CCDData
        Image array
    shift : donuts.image.Image
        The shift between the reference and this image is found
        in the image object for the shift measurment
    x : array
        X postions to extract photometry
    y : array
        Y postions to extract photometry
    rsi : array
        Inner sky annulii radii
    rso : array
        Outer sky annulii radii
    aper_radii : array
        Sequence of aperture radii to extract
    filename : str
        The name of the file we are extracting photometry from
    jd : astropy.time.Time
        JD at mid exposure
    bjd : astropy.time.Time
        BJD_TDB at mid exposure
    hjd : astropy.time.Time
        HJD as mid exposure
    phot_filename_prefix : str
        Prefix for the photometry output filename
        Default = 'rtp'
    ds9_name : str, optional
        Name of instance of the DS9 program
        Default = None
    draw_regions : bool, optional
        Boolean to set whether we draw the regions in DS9 or not
        Default = False
    gain : float, optional
        CCD gain setting
        Default = 1.00
    index_offset : int
        Index offset between python and ds9
        default = 1

    Returns
    -------
    None

    Raises
    ------
    None
    """
    x = np.array(x)-shift.x.value
    y = np.array(y)-shift.y.value
    # not CCDData object
    data = data.data
    # loop over each aperture size
    for r, r_aper in enumerate(aper_radii):
        # preamble
        out_str = "{0:s}  {1:.8f}  {2:.8f}  {3:.8f}  ".format(filename, jd.value, bjd.value, hjd.value)
        # do the phot
        flux, fluxerr, _ = sep.sum_circle(data, x, y, r_aper,
                                          subpix=0,
                                          bkgann=(rsi, rso),
                                          gain=gain)
        flux_w_sky, fluxerr_w_sky, _ = sep.sum_circle(data, x, y, r_aper,
                                                      subpix=0,
                                                      gain=gain)
        # loop over each object and compile the output
        for i, j, inner, outer, f, fe, fs, fse in zip(x, y, rsi, rso, flux, fluxerr, flux_w_sky, fluxerr_w_sky):
            if ds9_name and draw_regions and r == 0:
                annulus = '{{annulus {0} {1} {2} {3} # color=green}}'.format(i+index_offset,
                                                                             j+index_offset,
                                                                             inner, outer)
                circle = '{{circle {0} {1} {2} # colo=red}}'.format(i+index_offset,
                                                                    j+index_offset,
                                                                    r_aper)
                dset(ds9_name, f"regions command '{annulus}'")
                dset(ds9_name, f"regions command '{circle}'")
            max_pixel_value = find_max_pixel_value(data, int(i), int(j), r_aper+1)
            temp = "{:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  ".format(
                float(i), float(j), float(f), float(fe), float(fs), float(fse), float(max_pixel_value))
            out_str = out_str + temp
        out_str = out_str + "\n"

        # output the photometry with .photAPERTURE_SIZE extension
        with open(f"{phot_filename_prefix}.phot{r_aper}", 'a') as outfile:
            outfile.write(out_str)
