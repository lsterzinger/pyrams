@xr.register_dataset_accessor("RAMS")
class RAMSAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._lwc = None
        
    @property
    def lwc(self):
        self._lwc = (self._obj.RCP + self._obj.RDP + self._obj.RRP)
        return self._lwc
    
    @property
    def lwp(self):
        if self._lwc is None:
            self.lwc()
        self._lwp = self._lwc.integrate('z')
        return self._lwp