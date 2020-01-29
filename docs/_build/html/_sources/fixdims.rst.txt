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

RAMSlibs has two ways of overcoming this problem. The first is to rebuild the ``xr.Dataset`` with the correct dimenions 
using `fix_duplicate_dims <apiref.html#ramslibs.data_tools.fix_duplicate_dims>`_