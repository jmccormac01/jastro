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

### Example nightly config

### Example region file



# Utilities


# Contributors

James McCormac
