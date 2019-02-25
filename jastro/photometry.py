"""
Functions for extracting photometry
"""
import numpy as np
import sep
from jastro.ds9 import set_ds9
from jastro.reduction import find_max_pixel_value

# pylint: disable=invalid-name
# pylint: disable=no-member

def phot(data, shift, config, filename, jd, bjd, hjd,
         phot_filename_prefix="rtp", ds9_name=None, draw_regions=False,
         gain=1.00, index_offset=1.0):
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
    config : array-like
        Dictionary of the processing configuration
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

    Returns
    -------
    None

    Raises
    ------
    None
    """
    x = np.array(config['x'])-shift.x.value
    y = np.array(config['y'])-shift.y.value
    # not CCDData object
    data = data.data
    # output the photometry with .photAPERTURE_SIZE extension
    for r, r_aper in enumerate(config['aperture_r']):
        outfile = open("{0:s}.phot{1:d}".format(phot_filename_prefix, r_aper), 'a')
        # loop over each object in the file
        out_str = "{0:s}  {1:.8f}  {2:.8f}  {3:.8f}  ".format(filename, jd.value, bjd.value, hjd.value)
        for i, j, inner, outer in zip(x, y, config['inner'], config['outer']):
            if ds9_name and draw_regions and r == 0:
                annulus = '{{annulus {0} {1} {2} {3} # color=green}}'.format(i+index_offset,
                                                                             j+index_offset,
                                                                             inner, outer)
                circle = '{{circle {0} {1} {2} # colo=red}}'.format(i+index_offset,
                                                                    j+index_offset,
                                                                    r_aper)
                set_ds9(ds9_name, "regions command '{}'".format(annulus))
                set_ds9(ds9_name, "regions command '{}'".format(circle))
            flux, fluxerr, _ = sep.sum_circle(data, i, j, r_aper,
                                              subpix=0,
                                              bkgann=(inner, outer),
                                              gain=gain)
            flux_w_sky, fluxerr_w_sky, _ = sep.sum_circle(data, i, j, r_aper,
                                                          subpix=0,
                                                          gain=gain)
            max_pixel_value = find_max_pixel_value(data, int(i), int(j), r_aper+1)
            temp = "{:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  {:.2f}  ".format(
                float(i), float(j), float(flux), float(fluxerr),
                float(flux_w_sky), float(fluxerr_w_sky), float(max_pixel_value))
            out_str = out_str + temp
        out_str = out_str + "\n"
        outfile.write(out_str)
        outfile.close()
