# Astronomy code

# Installation

```sh
git clone https://github.com/jmccormac01/jastro.git
```

# Data Reduction and Aperture Photometry

The ```example.py``` script shows the steps for a typical
CCD reduction with aperture photometry to extract the
differential light curve of a varying source.

Two configuration files are needed:
   1. An instrument level config file (static parameters, see example below)
   1. A nightly config file with (variable parameters, see example below)

The config info is used to set up the data reduction and photometry routines.
For aperture photometry a region file corresponding to the reference image
must be supplied. The regions should mark the location of the sky annulii
for each comparison star (numbered 1, 2. ... N) and the target (mark with
a name, e.g. WASP-39). See example region file below

### Example instrument config

```json
{
    "filter": "RGB-R",
    "reference_image": "TOIXXX-CAM3-20230101200000.fits",
    "region_file": "TOIXXX-CAM3-20230101200000.reg",
    "binning": 10,
    "defocused": false,
    "ds9": true,
    "max_sep_shift": 2,
    "max_donuts_shift": 50,
    "output_reduced": false,
    "aperture_radii": [
        5.0,
        5.5,
        6.0,
        6.5,
        7.0
    ],
    "fit_parameters": [
        0.4,
        0.5
    ],
    "master_bias_filename": "master-bias-CAM3.fits",
    "master_dark_filename": "master-dark-CAM3.fits",
    "master_flat_filename": "master-flat-CAM3-RGB-R.fits"
}
```

### Example nightly config

```json
{
    "imager": {
        "border": 50,
        "filters": [
            "RGB-R",
            "RGB-G",
            "RGB-B"
        ],
        "ra_keyword": "MNTRA",
        "dec_keyword": "MNTDEC",
        "image_keyword": "SCIENCE",
        "flat_keyword": "FLAT",
        "bias_keyword": "BIAS",
        "dark_keyword": "DARK",
        "dateobs_start_keyword": "DATE-OBS",
        "exptime_keyword": "EXPTIME"
    },
    "observatory": {
        "olat": 28.760277,
        "olon": -17.879444,
        "elev": 2332
    },
    "sky": {
        "inner_radius": 20,
        "outer_radius": 25,
        "background_sigma": 3
    },
    "ds9": {
        "window_id": "LEOWASP",
        "index_offset": 1
    }
}
```

### Example region file

Contents of ```TOIXXX-CAM3-20230101200000.reg```

```
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
annulus(4812.2087,3195.6101,9.0,14) # text={TOIXXX}
annulus(4635.7325,3060.4728,20.0,30.0) # text={1}
annulus(4986.6286,3189.2174,20.0,30.0) # text={2}
annulus(5013.2188,3062.3943,20.0,30.0) # text={3}
annulus(4622.4847,3390.6434,20.0,30.0) # text={4}
annulus(5048.4832,3392.7397,20.0,30.0) # text={5}
```

# Utilities


# Contributors

James McCormac
