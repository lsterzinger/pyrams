"""Contains a collection of functions for derived variables."""
import numpy as np
from netCDF4 import Dataset as ncfile
from matplotlib import pyplot as plt
import xarray as xr
import warnings


class DataInfo():
    """
    A class to handle model data
    NOTE: Deprecated. Please use `DataVar()`
    
    ...

    Attributes
    ----------
    variable : str
        The name of the variable as found in the data files (e.g. "RTP")
    longname : str
        The long name of the variable (e.g. "Total Water Mixing Ratio")
    unit : str
        The unit of the variable (e.g. "kg/kg")
    data : ndarray
        The data array for the variable
    """

    def __init__(self, variable, longname, unit):
        warnings.warn('Note: data_tools.DataInfo is depricated.\
            Please move to data_tools.DataVar')
        self.variable = variable
        self.longname = longname
        self.unit = unit

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
        data : ndarray
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
    
    ...

    Parameters
    ----------
    varname :
        str The variable name as found in the data files (e.g. "RTP")
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
    data : ndarray
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
        for i in range(len(flist)):
            ds = xr.open_dataset(flist[i])
            self.data[i, :] = ds[self.varname]

    def purge_data(self):
        self.data = None


def flist_to_times(flist):
    """
    Creates a list of datetimes from a list of RAMS output
    variables.

    Note: This only works if there are no hyphens ('-') in the
    folder names.
.
    Parameters
    ----------
    flilst: list
        A list of files

    Returns
    -------
    times: list
        A list of times in datetime format.
    """
    from datetime import datetime

    times = np.zeros(len(flist), dtype=datetime)
    for i in range(len(flist)):
        traw = flist[i].split("-")[2:6]
        traw = "".join(traw)
        times[i] = datetime.strptime(traw, "%Y%m%d%H%M%S")

    return times


def vert_int(data, density, zheights, no_time=False):
    """
    Calculates the vertical integration given 3-D data, air density,
    and grid height

    Parameters
    -----------
    data: ndarray
        The data array. Can be in form (z, y, x), (z, x),
        (t, z, y, x), and (t, z, x)

    density: ndarray
        The air density, in (z, y, x) or (z, x)

    zheights: ndarray
        The height of the gridboxes, in (z, y, x) or (z, x)

    no_time: Bool, optional
        Flag for indicating lack of time dimension. Default=False

    Returns
    -------
    data_int: ndarray
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
    pressure : ndarray
         The 3-D pressure field (assumes time dimension, turn off
         with `no_time=True`)

    heights : ndarray
        The 3-D array of gridbox heights

    plevels : list
        List of pressure levels to interpolate to

    no_time=False: boolean
        Optional, set to `True` to indicate lack of time dimension.

    Returns
    -------
    press_height : ndarray
        The geopotential heights at the specified pressure levels
    """
    from metpy.interpolate import log_interpolate_1d

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
    topt : ndarray
        The 2-D topographic height information
        
    ztn : ndarray
        The ztn variable from the *head.txt files output from RAMS

    Returns
    -------
    z : ndarray
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


# def plot_domain2(variable, lats, lons, tmax):
#     from mpl_toolkits.basemap import Basemap

#     centlat = 38.2
#     centlon = -122.1
#     width = 925000
#     height = 700000
#     plt.ioff()
#     for t in range(0, tmax):
#         plt.figure(figsize=(12, 12))
#         m = Basemap(projection='stere', lon_0=centlon, lat_0=centlat,
#                     lat_ts=centlat, width=width, height=height)

#         m.drawcoastlines()
#         m.drawstates()
#         m.drawcountries()
#         parallels = np.arange(0., 90, 10.)
#         m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
#         meridians = np.arange(180., 360, 10.)
#         m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)

#         x, y = m(lons, lats)
#         m.contourf(x, y, variable[t, :, :])
#         plt.savefig(str(t) + ".png")


# def vert_int(variable, dimensions):
#     xmax = dimensions[0]
#     ymax = dimensions[1]
#     zmax = dimensions[2]
#     tmax = dimensions[3]

#     var_out = np.zeros((tmax, ymax, xmax))

#     for t in range(0, tmax):
#         for x in range(0, xmax):
#             for y in range(0, ymax):
#                 col_tot = 0

#                 for z in range(0, zmax):
#                     col_tot = col_tot + variable[t, z, y, x]
#                 var_out[t, y, x] = col_tot

#     return var_out


def z_levels_3d(ztn, topt):
    """
    Calculates the gridbox heights for a 3-D grid

    Parameters
    ----------
    ztn : list
        List of ztn values from RAMS *head.txt output

    topt : ndarray
        2-D array of topography height values

    Returns
    -------
    zheight : ndarray
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

    topt : ndarray
        1-D array of topography height values

    Returns
    -------
    zheight : ndarray
        2-D array of gridbox heights
    """

    xlen = topt.shape[1]
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
    pi : ndarray
        PI (modified exner function) from RAMS model output

    Returns
    -------
    pressure : ndarray
        Pressure (millibars)
    """

    p0 = 1000.
    cp = 1004.
    R = 287.

    return p0*np.power((pi/cp), (cp/R))


def temperature(theta, pi, degc=False):
    """
    Calculate the temperature from Exner function using THETA and PI.
    
    Parameters
    ----------
    theta : ndarray
        THETA (potential temperature) variable from RAMS output

    pi : ndarray
        PI (modified exner function) from RAMS output

    degc=False : boolean
        Optional, set to `True` to output temperature in Celsius instead of Kelvin

    Returns
    -------
    temperature : ndarray
        The temperature in Kelvin (or Celsius if `degc=True`)
    """
    cp = 1004.

    temp = theta*(pi/cp)
    if degc is True:
        temp = temp - 273.15

    return(temp)


def mslp(temp, press, height):
    """
    Calculate the mean sea level pressure

    Parameters
    ------
    temp : ndarray
        Temperature in Celsius

    press: ndarray
        Pressure in hPa

    height: ndarray
        height of the terrain in meters

    Returns
    -------
    p0 : ndarray
        The MSLP pressure (hPa)
    """
    # p0 = press*(1-(0.0065*height)/(temp + 0.0065*height + 273.15))**-5.257
    p0 = press * (1-(0.0065*height)/(temp+0.0065*height + 273.15))**-5.257

    return(p0)
