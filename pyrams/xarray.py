from logging import warning
from tokenize import String
from warnings import WarningMessage
import xarray as xr
@xr.register_dataset_accessor("rams")
class RAMSAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._lwc = None
        self._iwc = None

    @property
    def lwc(self):
        """
        Calculate liquid water content
        """

        vars = [
            'RCP',
            'RRP',
            'RDP'
        ]

        vsum = 0
        for v in vars:
            try:
                vsum = vsum + self._obj[v]
            except KeyError:
                print(f"Warning, {v} not found in dataset - skipping in LWC calculation")

        return (vsum) * self._obj.DN0
    
    @property
    def iwc(self):
        """
        Calculate Ice Water Content (kg/m3)
        """
        vars = [
            'RPP',
            'RSP',
            'RAP',
            'RGP',
            'RHP'
        ]
        vsum = 0
        for v in vars:
            try:
                vsum = vsum + self._obj[v]
            except KeyError:
                print(f"Warning, {v} not found in dataset - skipping in IWC calculation")

        return (vsum) * self._obj.DN0
    
    @property
    def lwp(self):
        """Calculate Liquid Water Path"""
        return self.lwc.integrate('z')

    @property
    def iwp(self):
        """Calculate Ice Water Path"""
        return self.iwc.integrate('z')

    def apply_variable_metadata(self, pint=True):
        """
        Applies metadata (unit and long_name) to RAMS output variables. 

        Parameters
        ----------
        pint: bool
            If true (default) will return with pint quantities added to variables via pint-xarray
        """
        import json
        import pkgutil

        ds = self._obj
        ramsvars = json.loads(pkgutil.get_data(__name__, 'rams-vars.json'))
        
        nokey = []
        for v in ds.variables:
            if v in ramsvars.keys():
                ds[v].attrs['units'] = ramsvars[v]['unit']
                ds[v].attrs['long_name'] = ramsvars[v]['long_name']
            else: 
                nokey.append(v)

        if len(nokey) > 0:
            print(f"Warning: no metadata found for variables {nokey}")

        if pint:
            import cf_xarray.units
            import pint_xarray
            return ds.pint.quantify()

    def fix_dims(
        self,
        flist=None, 
        dx=None,
        dz=None,
        z=None,
        dims = {
            'phony_dim_0' : 'x',
            'phony_dim_1' : 'y',
            'phony_dim_2' : 'z'
        },
        dt = None,
        ):
        """
        Calls pyrams.data_tools.create_xr_metadata()

        Parameters
        ----------  
        flist: List of file paths, optional
            List of filepaths, used to add datetimes to time dimension

        dt: ``np.timedelta64``
            Change the ``time`` coordinate to be a timedelta of unit ``dt``, 
            ``flist`` must be specified.

        dims: dict, optional
            Dict of dims to rename. defaults to ::

                dims = {
                    'phony_dim_0' : 'x',
                    'phony_dim_1' : 'y',
                    'phony_dim_2' : 'z'
                }

        dx: float, optional
            dx to add values to ``(x,y)`` dimensions
        
        dz: float, optional
            dz to add values to ``z`` dimension

        z: list, optional
            List of explicit ``z`` values to add to dimension

        Returns
        -------
        ds: ``xr.Dataset()``

        """

        from pyrams.data_tools import flist_to_times, create_xr_metadata
        import numpy as np

        return create_xr_metadata(
            self._obj, flist = flist, dims = dims, dx = dx,
            dz = dz, z = z, dt = dt)
