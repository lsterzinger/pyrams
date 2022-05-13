from logging import warning
import xarray as xr
@xr.register_dataset_accessor("rams")
class RAMSAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._lwc = None
        self._iwc = None

    @property
    def lwc(self):
        return (self._obj.RCP + self._obj.RDP + self._obj.RRP) * self._obj.DN0
    
    @property
    def iwc(self):
        return (self._obj.RPP + self._obj.RSP + self._obj.RAP + self._obj.RGP + self._obj.RHP) * self._obj.DN0
    
    @property
    def lwp(self):
        return self.lwc.integrate('z')

    @property
    def iwp(self):
        return self.iwc.integrate('z')

    def apply_metadata(self):
        import json

        ds = self._obj

        with open('./rams-vars.json', 'r') as inf:
            ramsvars = json.loads(inf.read())
        
        nokey = []
        for v in ds.variables:
            if v in ramsvars.keys():
                ds[v].attrs['unit'] = ramsvars[v]['unit']
                ds[v].attrs['long_name'] = ramsvars[v]['long_name']
            else: 
                nokey.append(v)

        if len(nokey) > 0:
            raise Warning(f"No metadata found for variables {nokey}")


    def fix_dims(
        self,
        flist=None, 
        dx=None,
        dims = {
            'phony_dim_0' : 'x',
            'phony_dim_1' : 'y',
            'phony_dim_2' : 'z'
        },
        dz=None,
        z=None):

        from pyrams.data_tools import flist_to_times
        import numpy as np

        ds = self._obj.copy()
        
        if flist and ds['time']:
            ds['time'] = flist_to_times(flist)
        
        if dx:
            Warning('dx')
            ds['x'] = np.arange(0, len(ds.x)) * dx
            ds['y'] = np.arange(0, len(ds.y)) * dx

        if dz:
            Warning('dz')
            ds['z'] = np.arange(0, len(ds.z))*dz 
        elif z:
            Warning('z')
            ds['z'] = z 

        return ds