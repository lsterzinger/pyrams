.. PyRAMS documentation master file, created by
   sphinx-quickstart on Wed Jan 29 09:11:31 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyRAMS Home
====================================
Welcome to the home of PyRAMS! This collection of tools was created by Lucas Sterzinger (lsterzinger@ucdavis.edu) to
ease the process of working with `RAMS <https://vandenheever.atmos.colostate.edu/vdhpage/rams.php>`_ model output. 
Check out the `code on GitHub! <https://github.com/lsterzinger/pyrams>`_


Installation
------------
PyRAMS can be installed  via conda::

   conda install -c conda-forge pyrams

or via pip::

   pip install pyrams

Interactive tutorial
--------------------
An interactive tutorial on PyRAMS integration with xarray can be found at this Binder link:

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/lsterzinger/pyrams/HEAD?labpath=examples%2Fxarray.ipynb

More interactive tutorials to follow (stay tuned!)

Basic usage
-----------
The biggest feature of PyRAMS is its extension of `xarray <https://docs.xarray.dev/en/stable/>`_. 
This functionality is outlined in :doc:`xarray_metadata` and an interactive example can be run 
in-browser via the Binder link above.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
    
   xarray_metadata
   fixdims
   api/index
