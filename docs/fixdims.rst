Fixing dimensions
=================
RAMS output, by default, uses dimensions named ``phony_dim_0``, ``phony_dim_1``, etc. The problem lies when two dimensions
have the same size. For example below, the data has the dimensions ``['z', 'y', 'x'] = [300, 32, 32]    ``. 
.. code-block:: 

    netcdf c-A-2008-08-31-000000-g1 {
    dimensions:
            phony_dim_0 = 32 ;
            phony_dim_1 = 200 ;
    variables:
            float AGGREGATET(phony_dim_1, phony_dim_0, phony_dim_0) ;
    }

However, since ``x=y=32``, RAMS assigns the dimension ``phony_dim_0`` to both x `and` y, causing errors notably in 
`xr.open_mfdataset() <http://xarray.pydata.org/en/stable/generated/xarray.open_mfdataset.html>`_

PyRAMS has two ways of overcoming this problem. The first is to rebuild the ``xr.Dataset`` with the correct dimenions 
using `fix_duplicate_dims <apiref.html#pyrams.data_tools.fix_duplicate_dims>`/


.. code-block:: python

    from pyrams.data_tools import fix_duplicate_dims

    ds = xr.open_dataset('./dataset.h5')

    # This line will replace 'phony_dim_0' with 'x' and 'y'
    ds_new = fix_duplicate_dims(ds, ['y', 'x'], 'phony_dim_0')


You can also rewrite a list of RAMS output files into netCDF with renamed dimensions:

.. code-block:: python

    from pyrams.data_tools import rewrite_to_netcdf
    from glob import glob

    flist = glob('/path/to/files/*.h5')
    rewrite_to_netcdf(flist, '/path/to/write/files/', ['y', 'x'], 'phony_dim_0', prefix='dimfix') 
    # prefix is optional, and defaults to 'dimfix'
