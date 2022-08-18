"""Contains a collections of functions for working with RAMS data in Python"""
from operator import concat
import numpy as np
from netCDF4 import Dataset as ncfile
from matplotlib import pyplot as plt
import xarray as xr
import warnings
from tqdm import tqdm
import pandas as pd
from datetime import datetime
from metpy.interpolate import log_interpolate_1d
import re

class DataInfo():
    """
    Deprecated. Please use `DataVar()`

    A class to handle model data

    Attributes
    ----------
    variable : str
        The name of the variable as found in the data files (e.g. "RTP")
    longname : str
        The long name of the variable (e.g. "Total Water Mixing Ratio")
    unit : str
        The unit of the variable (e.g. "kg/kg")
    data : numpy.ndarray
        The data array for the variable
    """

    def __init__(self, variable, longname, unit):
        warnings.warn('Note: data_tools.DataInfo is depricated.\
            Please move to data_tools.DataVar')
        self.variable = variable
        self.longname = longname
        self.unit = unit
        self.data = None

    def get_data(self, datadir, simulation):
        """
        Parameters
        ----------
        datadir : str
            The path of the data files
        simulation : str
            Name of the subfolder that the data is found in
            (e.g. "feb2014_control")

        Returns
        -------
        data : numpy.ndarray
            The data for the desired variable
        """

        file = xr.open_mfdataset(datadir + simulation + "*g2.h5",
                                 concat_dim='TIME')
        data = file[self.variable]
        return data


class DataVar():
    """
    A new class created to manage variables, their names, and units
    (replaces `DataInfo`)

    Parameters
    ----------
    varname : str
        The variable name as found in the data files (e.g. "RTP")
    longname : str
        Optional. The long name of the variable
        (e.g. "Total Water Mixing Ratio")
    unit : str
        Optional. The unit of the variable (e.g. "kg/kg")

    Attributes
    ----------
    varname :
        str The variable name as found in the data files (e.g. "RTP")
    longname : str
        The long name of the variable (e.g. "Total Water Mixing Ratio")
    unit : str
        The unit of the variable (e.g. "kg/kg")
    data : numpy.ndarray
        The data for the variable
    """

    def __init__(self, varname, longname=None, unit=None):
        self.varname = varname
        self.unit = unit
        self.longname = longname

    def get_data(self, flist):
        """
        Pulls data from a list of files (flist) and puts it into a single array

        Arguments
        ---------
        flist : list
            A list of sorted file paths
        """
        print(f'Opening {self.varname}')

        # Get dimensions
        dims = list(xr.open_dataset(flist[0])[self.varname].shape)
        dims.insert(0, len(flist))

        # Create empty array
        self.data = np.zeros(dims)

        # Get data
        for i in tqdm(range(len(flist))):
            ds = xr.open_dataset(flist[i])
            self.data[i, :] = ds[self.varname]

    def purge_data(self):
        self.data = None


def domain_mean_netcdf(ds_with_metadata, outfile, vars=None):
    """
    Writes x/y domain-average from from an xarray dataset to `outfile` as NetCDF.

    Arguments
    ---------
    ds_with_metadata : xr.Dataset
        An xarray dataset created with `pyrams.datatools.create_xr_dataset()`

    outfile : str
        Name of output file
    
    vars : list (optional)
        List of variable names to write. Default is to process and write all variables
    """

    from os.path import exists
    from tqdm import tqdm

    ds = ds_with_metadata

    if vars is None:
        vars = [i for i in ds.data_vars]

    if exists(outfile):
        raise Exception(f"Error: {outfile} already exists")
                        
    for v in tqdm(vars):
        try:

            a = ds[v].mean(dim=('x', 'y'))
            if exists(outfile):
                mode = 'a'
            else:
                mode = 'w'
            a.to_netcdf(outfile, mode=mode)
            # print(a)
        except ValueError:
            print(f"Variable {f} ")
            continue

def rewrite_to_netcdf(flist, output_path, duped_dims, phony_dim, prefix='dimfix', single_file=False, compression_level=None):
    """
    Rewrites RAMS standard output files as netCDF4 with fixed dimension data, using 
    ``data_tools.fix_duplicate_dims()``

    Arguments
    ---------
    flist : list of str
        A list of file paths (recommend using ``sorted(glob.glob('/path/to/files/*g1.h5'))`` or similar)

    output_path : str
        Path where new files will be written

    duped_dims : list of str
        List of dimensions that are duplicated, in order (e.g. ``['y', 'x']``)

    phony_dim : string
        Name of duplicate dimension in `ds`, often ``phony_dim_0``

    prefix : string
        Prefix for output files, defaults to `dimfix`

    single_file : bool (optional)
        If `True`, will combine all files into a single file with name `<prefix>.nc`. Defaults to `False`.

    compression_level : int
        If specified, data will be compressed on the given level (0-9 are valid). Defaults to `None`
    """
    dslist = []
    for f in tqdm(flist):
        str_date = '-'.join(f.split('/')[-1].split("-")[2:6])
        date = np.datetime64(pd.to_datetime(str_date))

        if duped_dims:
            ds = fix_duplicate_dims(xr.open_dataset(f), duped_dims, phony_dim)

        else:
            ds = xr.open_dataset(f)

        ds = ds.expand_dims({'time': [date]})

        if not single_file:
            ds.to_netcdf(f'{output_path}/{prefix}_{str_date}.nc',
                         unlimited_dims=['time'])

        if single_file:
            dslist.append(ds)

        ds.close()

    if single_file:
        ds = xr.concat(dslist, dim='time')
        encoding = None

        # If compression_level is defined, compress files
        if compression_level:

            # Catch invalid values
            if type(compression_level) is not int or compression_level < 0 or compression_level > 9:
                raise ValueError(
                    "compression_level must be an integer between 0 and 9 (inclusive)")

            comp = dict(zlib=True, complevel=compression_level)
            encoding = {var: comp for var in ds.data_vars}

        ds.to_netcdf(f"{output_path}/{prefix}.nc", encoding=encoding)
    return


def fix_duplicate_dims(ds, duped_dims, phony_dim):
    """
    Fixes duplicate dimensions (often with the same amount of `x` and `y` gridpoints), 
    for use with xarray.open_mfdataset and xarray.combine_nested.

    Arguments
    ---------
    ds : xarray.Dataset
        The dataset to be fixed

    duped_dims : list of str
        List of dimensions that are duplicated, in order (e.g. `['y', 'x']`)

    phony_dim : string
        Name of duplicate dimension in `ds`, often `'phony_dim_0'`


    Returns
    -------
    ds_new : xarray.Dataset
        New dataset with fixed dimension names
    """
    dims = dict(ds.dims)

    try:
        dupe_dim = dims[phony_dim]
    except KeyError:
        print(f'Error, duplicate dimension must be \'phony_dim_0\'')
        return

    dims.pop(phony_dim)

    for d in duped_dims:
        dims[d] = dupe_dim

    ds_new = xr.Dataset()

    for v in ds.variables:
        dvar = ds[v]

        # Check if phony_dim exists in variable dimensions
        if phony_dim in dvar.dims:
            vardims = list(dvar.dims)
            indices = [i for i, x in enumerate(vardims) if x == phony_dim]
            for i, ind in enumerate(indices):
                vardims[ind] = duped_dims[i]
            vardims = tuple(vardims)

            ds_new[v] = (vardims, ds[v])

        else:
            vardims = dvar.dims
            ds_new[v] = (vardims, ds[v])

    ds.close()
    return(ds_new)


def flist_to_times(flist):
    """
    Creates a list of datetimes from a list of RAMS output
    variables.

    Function uses regex to find the pattern "YYYY-mm-dd-HHMMSS" in
    the file path and converts to a np.datetime64 object.

    Parameters
    ----------
    flilst: list
        A list of files

    Returns
    -------
    times: list
        A list of times in np.datetime64 format.
    """
    dtregex = r"[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{6}"

    times = np.zeros(len(flist), dtype=datetime)
    for i,f in enumerate(flist):

        tarr = re.findall(dtregex, f)
        if tarr == []:
            raise SyntaxError(f"No datetimes of form \"YYYY-MM-DD-HHMMSS\" were found in {f}")

        if len(tarr) == 1:
            traw = tarr[0]
        else:
            raise SyntaxError(f"More than one date found in path {f}")

        traw = tarr[0]
        times[i] = datetime.strptime(traw, "%Y-%m-%d-%H%M%S")

    return times.astype(np.datetime64)


def create_xr_metadata(
    ds,
    flist = None,
    dims = {
        'phony_dim_0' : 'x',
        'phony_dim_1' : 'y',
        'phony_dim_2' : 'z'
    },
    dx = None,
    dz = None,
    z = None,
    dt = None,
):
    import numpy as np
    """
    Adds metadata to ``xr.Dataset()``.

    Parameters
    ----------
    ds: ``xarray.Dataset``
        Dataset
    
    flist: List of file paths, optional
        List of filepaths, used to add datetimes to time dimension

    dt: ``String``
        One of ``['second', 'minute', 'hour', 'day']``
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

    
    try:
        ds = ds.rename(dims)
    except ValueError:
        print("Dimensions have already been renamed - skipping and continuing with other metadata")

    if flist:
        if type(flist) is str:
            flist = [flist]
        
        try:
            ds['time'] = flist_to_times(flist)
        except KeyError:
            ds = ds.assign(time=flist_to_times(flist))

    dt_options = {
        'second' : np.timedelta64(1, 's'),
        'minute' : np.timedelta64(1, 'm'),
        'hour' : np.timedelta64(1, 'h'),
        'day' : np.timedelta64(1, 'D'),
    }
    if dt in dt_options:

        # If we're getting timestamps from `flist`, convert to timedelta 
        # and divide by dt to get in specified units
        if flist:
            t = (ds['time'] - ds['time'][0]) / dt_options[dt]
            ds['time'] = t
        
        else:
            raise TypeError('`flist` must be definied to use `dt`.')
        # If not, assume dt describes length between time indices
        # else:
        #     ds['time'] = np.arange(0, dt*len(ds['time']), dt)
    

        # Add unit
        ds['time'].attrs['unit'] = dt

    # If not in dt_options or None, raise error
    elif dt is not None:
        raise TypeError(f"`dt` must one of {list(dt_options.keys())}")

    if dx:
        ds['x'] = np.arange(0, len(ds.x)) * dx
        ds['y'] = np.arange(0, len(ds.y)) * dx

        ds['x'].attrs = {'units' : 'm'}
        ds['y'].attrs = {'units' : 'm'}
    
    if dz:
        ds['z'] = np.arange(0, len(ds.z)) * dz
        ds['z'].attrs = {'units' : 'm'}

    if z:
        ds['z'] = z
        ds['z'].attrs = {'units' : 'm'}

    return ds


def habit_count(habits, tmax):
    """
    Takes 3D habit data and tmax (number of time steps) and returns the number
    of each habit at each time step.
    """
    # Initialize empty array
    count = np.zeros((tmax, 12))

    for time in range(0, tmax):
        unique, counts = np.unique(habits[time, :, :, :],
                                   return_counts=True)

        for i in range(0, len(unique)):
            if np.ma.is_masked(unique[i]) is False:
                count[time, int(unique[i])] = counts[i]

    return count


def press_level(pressure, heights, plevels, no_time=False):
    """
    Calculates geopotential heights at a given pressure level

    Parameters
    ----------
    pressure : numpy.ndarray
         The 3-D pressure field (assumes time dimension, turn off
         with `no_time=True`)

    heights : numpy.ndarray
        The 3-D array of gridbox heights

    plevels : list
        List of pressure levels to interpolate to

    no_time=False: bool 
        Optional, set to `True` to indicate lack of time dimension.

    Returns
    -------
    press_height : numpy.ndarray
        The geopotential heights at the specified pressure levels
    """

    if no_time is False:
        try:
            tlen, zlen, ylen, xlen = pressure.shape

            press_height = np.zeros((tlen, ylen, xlen))
            for t in range(0, tlen):
                for x in range(0, xlen):
                    for y in range(0, ylen):
                        press_height[t, y, x] =\
                            log_interpolate_1d(plevels, pressure[t, :, y, x],
                                               heights[:, y, x])
        except ValueError:
            print("Error in dimensions, trying with no_time=True")
            no_time = True

    elif no_time is True:
        try:
            xlen, ylen, xlen = pressure.shape

            press_height = np.zeros((ylen, xlen))
            for x in range(0, xlen):
                for y in range(0, ylen):
                    press_height[t, y, x] =\
                        log_interpolate_1d(plevels, pressure[t, :, y, x],
                                           heights[:, y, x])
        except ValueError:
            print("Error in dimensions")

    return press_height


def calc_height(topt, ztn):
    """
    Calculates the height of each grid box

    Parameters
    ----------
    topt : numpy.ndarray
        The 2-D topographic height information

    ztn : numpy.ndarray
        The ztn variable from the *head.txt files output from RAMS

    Returns
    -------
    z : numpy.ndarray
        A 3-D array of the heights of each gridbox
    """
    ylen, xlen = topt.shape
    zlen = len(ztn)

    z = np.zeros((zlen, ylen, xlen))
    ztop = ztn[59]
    for x in range(0, xlen):
        for y in range(0, ylen):
            z[:, y, x] = ztn * (1 - (topt[y, x]/ztop)) + topt[y, x]
    return z


def z_levels_3d(ztn, topt):
    """
    Calculates the gridbox heights for a 3-D grid

    Parameters
    ----------
    ztn : list
        List of ztn values from RAMS *head.txt output

    topt : numpy.ndarray
        2-D array of topography height values

    Returns
    -------
    zheight : numpy.ndarray
        3-D array of gridbox heights
    """

    ylen, xlen = topt.shape
    zlen = len(ztn)

    zheight = np.zeros((zlen, ylen, xlen))
    for x in range(0, xlen):
        for y in range(0, ylen):
            for z in range(0, zlen):
                zheight[z, y, x] = ztn[z] * \
                    (1-topt[y, x]/ztn[zlen-1])+topt[y, x]

    return zheight


def z_levels_2d(ztn, topt):
    """
    Calculates the gridbox heights for a 2-D grid

    Parameters
    ----------
    ztn : list
        List of ztn values from RAMS *head.txt output

    topt : numpy.ndarray
        1-D array of topography height values

    Returns
    -------
    zheight : numpy.ndarray
        2-D array of gridbox heights
    """

    xlen = topt.shape[0]
    zlen = len(ztn)

    zheight = np.zeros((zlen, xlen))
    for x in range(0, xlen):
        for z in range(0, zlen):
            zheight[z, x] = ztn[z] * \
                (1-topt[x]/ztn[zlen-1])+topt[x]

    return zheight


def build_mfdataset(path, **kwargs):
    """
    Build and xarray dataset with a time dimension

    Parameters
    ----------
    path : string
        Path to folder containing files

    **kwargs
        Additional arguments to pass to xarray
    """
    from glob import glob
    flist = sorted(glob(path + '/*.h5'))
    print(f'Building dataset with {len(flist)} files')

    return xr.open_mfdataset(flist, combine='nested', concat_dim='time', **kwargs)