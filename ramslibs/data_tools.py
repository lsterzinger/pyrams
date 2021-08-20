"""Contains a collections of functions for working with RAMS data in Python"""
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
    z = None
):
    """
    Adds metadata to ``xr.Dataset()``.

    Parameters
    ----------
    ds: ``xarray.Dataset``
        Dataset
    
    flist: List of file paths, optional
        List of filepaths, used to add datetimes to time dimension

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
    
    ds = ds.rename(dims)

    if flist is not None:
        ds['time'] = flist_to_times(flist)

    if dx is not None:
        ds['x'] = np.arange(0, len(ds.x)) * dx
        ds['y'] = np.arange(0, len(ds.y)) * dx
    
    if dz is not None:
        ds['z'] = np.arange(0, len(ds.z)) * dz

    if z is not None:
        ds['z'] = x

    return ds



def vert_int(data, density, zheights, no_time=False):
    """
    Calculates the vertical integration given 3-D data, air density,
    and grid height

    Parameters
    -----------
    data: numpy.ndarray
        The data array. Can be in form (z, y, x), (z, x),
        (t, z, y, x), and (t, z, x)

    density: numpy.ndarray
        The air density, in (z, y, x) or (z, x)

    zheights: numpy.ndarray
        The height of the gridboxes, in (z, y, x) or (z, x)

    no_time: bool, optional
        Flag for indicating lack of time dimension. Default=False

    Returns
    -------
    data_int: numpy.ndarray
        The vertically integrated array. Same dimensions as `data` except
        without the 'z' dimension
    """

    # if 3-D, assume z,y,x integrate to y,x
    # If there's no time coordinate:
    if no_time:
        if len(data.shape) == 3:  # (z, y, x)
            zax, yax, xax = 0, 1, 2
            zlen, ylen, xlen = list(data.shape)
            data_int = np.zeros((ylen, xlen))

        elif len(data.shape) == 2:  # (z, x)
            zax, xax = 0, 1
            zlen, xlen = list(data.shape)
            data_int = np.zeros((xlen))

        else:
            raise ValueError("Error with dimensions...\
                are you sure there's no time dimension?")

        # Do vertical integration
        data = data * density  # Multiply by density to get kg/m^3
        for z in range(zlen-1):
            data_int += data[z, :] * (zheights[z+1, :] - zheights[z, :])

    # If there is a time coordinate
    else:
        if len(data.shape) == 4:  # (t, z, y, x)
            tax, zax, yax, xax = 0, 1, 2, 3
            tlen, zlen, ylen, xlen = list(data.shape)
            data_int = np.zeros((tlen, ylen, xlen))

        elif len(data.shape) == 3:  # (t, z, x)
            tax, zax, xax = 0, 1, 2
            tlen, zlen, xlen = list(data.shape)
            data_int = np.zeros((tlen, xlen))

        else:
            raise ValueError("Error with dimensions,\
                are you sure you mean to have a time dimension?")

        # Do vertical integration
        for t in range(tlen):
            data[t, :] = data[t, :] * density  # Multiply by density to kg/m^3
            for z in range(zlen-1):
                data_int[t, :] += data[t, z, :] *\
                    (zheights[z+1, :] - zheights[z, :])

    return data_int


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


def pressure(pi):
    """
    Calculate the pressure from Exner function using PI.

    Parameters
    ----------
    pi : numpy.ndarray
        PI (modified exner function) from RAMS model output

    Returns
    -------
    pressure : numpy.ndarray
        Pressure (millibars)
    """

    p0 = 1000.
    cp = 1004.
    R = 287.

    return p0*np.power((pi/cp), (cp/R))
