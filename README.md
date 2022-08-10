# PyRAMS - Useful tools for working with RAMS data

[![Documentation Status](https://readthedocs.org/projects/pyrams/badge/?version=latest)](https://pyrams.readthedocs.io/en/latest/?badge=latest)[![DOI](https://zenodo.org/badge/176599749.svg)](https://zenodo.org/badge/latestdoi/176599749)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyrams/badges/version.svg)](https://anaconda.org/conda-forge/pyrams)

## Basic Usage
The biggest feature of PyRAMS is its extension of [xarray](https://docs.xarray.dev/en/stable/). 
PyRAMS adds named dimensions, coordinate variables, and variable metadata to xarray datasets created from RAMS data.
This functionality is [outlined here](https://pyrams.readthedocs.io/en/latest/xarray_metadata.html) and an interactive example can be [run 
in-browser here](https://mybinder.org/v2/gh/lsterzinger/pyrams/HEAD?labpath=examples%2Fxarray.ipynb). Other functions can be found in the [API reference](https://pyrams.readthedocs.io/en/latest/api/index.html).

***
Example notebook to show how to use PyRAMS with xarray can be run inside the browser here:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lsterzinger/pyrams/HEAD?labpath=examples%2Fxarray.ipynb)

***
Documentation can be [found here](https://pyrams.readthedocs.io/en/latest)

Install from conda-forge:
```
conda install -c conda-forge pyrams
```

or with pip:
```
pip install pyrams
``` 
PyRAMS is distributed under an [MIT License](LICENSE)
