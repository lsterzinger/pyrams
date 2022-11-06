import xarray as xr
import logging


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
                print(
                    f"Warning, {v} not found in dataset - skipping in LWC calculation")

        lwcout = (vsum) * self._obj.DN0
        lwcout.attrs['long_name'] = 'Liquid Water Content'
        return lwcout

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
                print(
                    f"Warning, {v} not found in dataset - skipping in IWC calculation")

        iwcout = (vsum) * self._obj.DN0
        iwcout.attrs['long_name'] = 'Ice Water Content'
        return iwcout

    @property
    def lwp(self):
        """Calculate Liquid Water Path"""
        logging.debug('LWP calculation')
        lwpout = self.lwc.integrate('z')
        lwpout.attrs['long_name'] = 'Liquid Water Path'

        # pint-xarray does not consider units of coordinates, so if we're tracking units
        # then multiply by meters to adjust for integration
        try:
            import pint
            logging.debug(f'LWP units: {lwpout.pint.units}')
            if lwpout.pint.units == pint.Unit('kg/m^3'):
                logging.debug('Units in LWP detected, manually adjusting')
                lwpout = lwpout.pint.dequantify()
                lwpout = lwpout.pint.quantify('kilogram / meter ** 2')
            return lwpout
        except AttributeError as a:
            print(a)
            logging.debug("No LWP units detected")
            return lwpout

    @property
    def iwp(self):
        """Calculate Ice Water Path"""
        logging.debug('IWP calculation')
        iwpout = self.iwc.integrate('z')
        iwpout.attrs['long_name'] = 'Ice Water Path'

        # pint-xarray does not consider units of coordinates, so if we're tracking units
        # then multiply by meters to adjust for integration
        try:
            import pint
            logging.debug(f'IWP units: {iwpout.pint.units}')
            if iwpout.pint.units == pint.Unit('kg/m^3'):
                logging.debug('Units in IWP detected, manually adjusting')
                iwpout = iwpout.pint.dequantify()
                iwpout = iwpout.pint.quantify('kilogram / meter ** 2')
            return iwpout
        except AttributeError as a:
            print(a)
            logging.debug("No IWP units detected")
            return iwpout

    @property
    def cloudradius(self):
        """Calculate cloud droplet mean size"""
        from .thermo import rho_w
        import numpy as np

        ds = self._obj
        rcp = ds.RCP
        ccp = ds.CCP

        r = ((3/4) * ((rcp/ccp)/(np.pi * rho_w)))**(1/3)
        r.attrs['long_name'] = 'Mean Cloud Droplet Radius'
        r.attrs['units'] = 'm'
        try:
            r = r.pint.dequantify()
            r = r.pint.quantify('m')
        except AttributeError:
            print("No Pint units detected, output will be in meters")
            r.attrs['units'] = 'm'
        return r

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
        dims={
            'phony_dim_0': 'x',
            'phony_dim_1': 'y',
            'phony_dim_2': 'z'
        },
        dt=None,
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
            self._obj, flist=flist, dims=dims, dx=dx,
            dz=dz, z=z, dt=dt)
