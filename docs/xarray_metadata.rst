Xarray Metadata
===============

RAMS outputs in hdf5 format, which does not have named dimensions. ``RAMSLibs.data_tools.create_xr_metadata()`` created named dimensions 
and populates them with the correct values. For example:

.. code-block:: python

    import xarray as xr
    from ramslibs.data_tools import create_xr_metadata
    from glob import glob

    # Get a list of files
    flist = sorted(glob('/home/lsterzin/arctic_model/final_sims/oliktok/output/control_newnudge/*.h5'))

    # Open the first 10 files for this example
    ds = xr.open_mfdataset(flist[:10], combine='nested', concat_dim='time')

    ds = create_xr_metadata(
        ds, 
        flist=flist[:10],
        dx = 62.5,
        dz = 6.26
    )

The ``create_xr_metadata`` function:
    * Renames ``phony_dim_0``, ``phony_dim_1``, and ``phony_dim_2`` to ``x``, ``y``, and ``z``
    * Uses ``flist`` to generate datetimes from the filenames, and adds them to the ``time`` coordinate
    * Uses ``dx`` to generate x/y coordinates
    * Uses ``dz`` to generate z coordinates
        - If using variably-spaced vertical coordinates, you can also pass a list of ``z`` values with ``z = [0, 10, 20, ...]``