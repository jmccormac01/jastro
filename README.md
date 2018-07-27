# Astronomy code

# Installation

```sh
git clone https://github.com/jmccormac01/jastro.git
```

# Astronomical constants

A list of useful constants which appear frequently in astronomical calculations

### Accessing constants

```python
>>> from jastro.constants import Constants as c
>>> print(c.G)
6.67408e-11
```

### Listing supported constants

```python
>>> c._list_constants()
G Msun Rsun Lsun Mjup Rjup Mproton boltzmann AU sec_d
```

# Light curve handling

A series of scripts for processing photometric light curves

### Converting differential magnitudes to relative fluxes



# Contributors

James McCormac
