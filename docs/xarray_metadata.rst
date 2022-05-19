Xarray Compatibility
====================

PyRAMS includes an xarray accessor to apply units and add other metadata to xarray datasets
built on RAMS data.

For example - variable metadata can be added with:

.. code-block:: python

    import xarray as xr
    from pyrams.xarray
    from glob import glob

    # Get a list of files
    flist = sorted(glob('/path/to/data/*.h5'))

    ds = xr.open_mfdataset(flist, combine='nested', concat_dim='time')

    ds.rams.apply_variable_metadata()

The ``ds.rams.apply_variable_metadata()`` function adds units and long names to variable attributes in xarray. 
The variable metadata was compiled from the `RAMS Variable List <https://vandenheever.atmos.colostate.edu/vdhpage/rams/docs/RAMS-VariableList.pdf>`_
and was converted to a `JSON file included within the PyRAMS source code <https://github.com/lsterzinger/pyrams/blob/main/pyrams/rams-vars.json>`_.

This can also be used to apply ``pyrams.data_tools.create_xr_metadata``:

.. code-block:: python

    ds = ds.rams.fix_dims(flist=flist, dx=62.5, dz=6.25)


This function
    * Renames ``phony_dim_0``, ``phony_dim_1``, and ``phony_dim_2`` to ``x``, ``y``, and ``z``
    * Uses ``flist`` to generate datetimes from the filenames, and adds them to the ``time`` coordinate
    * Uses ``dx`` to generate x/y coordinates
    * Uses ``dz`` to generate z coordinates
        - If using variably-spaced vertical coordinates, you can also pass a list of ``z`` values with ``z = [0, 10, 20, ...]``


Running both of the above tools on a dataset will yield more metadata when looking at variables, such as ``ds.RCP`` (cloud water content):

.. code-block:: parsed-literal

    <xarray.DataArray 'RCP' (z: 200, x: 96, y: 96)>
    [1843200 values with dtype=float32]
    Coordinates:
    * x        (x) float64 0.0 62.5 125.0 187.5 ... 5.812e+03 5.875e+03 5.938e+03
    * y        (y) float64 0.0 62.5 125.0 187.5 ... 5.812e+03 5.875e+03 5.938e+03
    * z        (z) float64 0.0 6.25 12.5 18.75 ... 1.231e+03 1.238e+03 1.244e+03
    Attributes:
        unit:       kg/kg
        long_name:  cloud mixing ratio