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