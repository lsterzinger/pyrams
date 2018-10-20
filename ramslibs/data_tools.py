r"""Contains a collection of functions for derived variables."""
import numpy as np
from netCDF4 import Dataset as ncfile
from matplotlib import pyplot as plt


def habit_count(habits, tmax):
    r"""
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


def press_level(pressure, heights, plevels, xyt_dimensions):
    r"""Returns geopotential heights at a given pressure level"""
    from metpy.interpolate import log_interpolate_1d
    from metpy.units import units

    xmax = xyt_dimensions[0]
    ymax = xyt_dimensions[1]
    tmax = xyt_dimensions[2]

    press_height = np.zeros((tmax, ymax, xmax)) * units.meter

    for t in range(0, tmax):
        for x in range(0, xmax):
            for y in range(0, ymax):
                press_height[t, y, x] =\
                 log_interpolate_1d(plevels, pressure[t, :, y, x],
                                    heights[:, y, x])
    return press_height


def calc_height(topt, dimensions):
    xmax = dimensions[0]
    ymax = dimensions[1]
    zmax = dimensions[2]

    z = np.zeros((zmax, ymax, xmax))
    ztop = units.ztn[59]
    for x in range(0, xmax):
        for y in range(0, ymax):
            z[:, y, x] = units.ztn * (1 - (topt[y, x]/ztop)) + topt[y, x]
        return z


def plot_domain2(variable, lats, lons, tmax):
    from mpl_toolkits.basemap import Basemap

    centlat = 38.2
    centlon = -122.1
    width = 925000
    height = 700000
    plt.ioff()
    for t in range(0, tmax):
        plt.figure(figsize=(12, 12))
        m = Basemap(projection='stere', lon_0=centlon, lat_0=centlat,
                    lat_ts=centlat, width=width, height=height)

        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
        parallels = np.arange(0., 90, 10.)
        m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
        meridians = np.arange(180., 360, 10.)
        m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)

        x, y = m(lons, lats)
        m.contourf(x, y, variable[t, :, :])
        plt.savefig(str(t) + ".png")


def vert_int(variable, dimensions):
    xmax = dimensions[0]
    ymax = dimensions[1]
    zmax = dimensions[2]
    tmax = dimensions[3]

    var_out = np.zeros((tmax, ymax, xmax))

    for t in range(0, tmax):
        for x in range(0, xmax):
            for y in range(0, ymax):
                col_tot = 0

                for z in range(0, zmax):
                    col_tot = col_tot + variable[t, z, y, x]
                var_out[t, y, x] = col_tot

    return var_out


def import_variable(rootpath, habit, varname, domain_number, dimensions):
    import glob
    filepath = rootpath + "feb2014_" + habit
    t = 0

    if(len(dimensions) == 3):
        xmax = dimensions[0]
        ymax = dimensions[1]
        tmax = dimensions[2]
        variable = np.zeros((tmax, ymax, xmax))

        files = []

        for name in glob.glob(filepath+'/*-g'+str(domain_number)+".h5"):
            if t == tmax:
                break
            files.append(name)
            t = t + 1

        files.sort()
        t = 0

        for name in files:
            input_file = ncfile(name, 'r')
            print name
            variable[t, :, :] = input_file[varname][:, :]
            t = t + 1
            input_file.close()

    elif(len(dimensions) == 4):
        xmax = dimensions[0]
        ymax = dimensions[1]
        zmax = dimensions[2]
        tmax = dimensions[3]
        variable = np.zeros((tmax, zmax, ymax, xmax))

        files = []

        for name in glob.glob(filepath+'/*-g'+str(domain_number)+".h5"):
            if t == tmax:
                break
            files.append(name)
            t = t+1

        files.sort()
        t = 0

        for name in files:
            input_file = ncfile(name)
            print name
            variable[t, :, :, :] = input_file[varname][:]
            t = t + 1
            input_file.close()

    return variable


def domain_average_2d(variable, dimensions):
    # Do domain average
    xmax = dimensions[0]
    ymax = dimensions[1]
    tmax = dimensions[2]
    average = np.zeros(tmax)

    for t in range(0, tmax):
        sum = 0
        for x in range(0, xmax):
            for y in range(0, ymax):
                sum = sum + variable[t, y, x]

        average[t] = sum/(xmax*ymax)
    return average


def time_average_3d(variable, tmax, zmax, ymax, xmax):
    total = 0
    output_var = np.zeros(tmax)

    for t in range(0, tmax):
        for x in range(0, xmax):
            for y in range(0, ymax):
                for z in range(0, zmax):
                    total = total + variable[t, z, y, x]
        output_var[t] = total/(xmax*ymax)*zmax
        total = 0
    return output_var


def time_average_2d(variable, tmax, ymax, xmax):
    total = 0
    output_var = np.zeros(tmax)

    for t in range(0, tmax):
        for x in range(0, xmax):
            for y in range(0, ymax):
                total = total + variable[t, y, x]
        output_var[t] = total/(xmax*ymax)
        total = 0
    return output_var


def z_levels_3d(ztn, topt, xmax, ymax, zmax):
    zheight = np.zeros((zmax, ymax, xmax))
    for x in range(0, xmax):
        for y in range(0, ymax):
            for z in range(0, zmax):
                zheight[z, y, x] = ztn[z] * \
                    (1-topt[y, x]/ztn[zmax-1])+topt[y, x]

    return zheight


def z_levels_2d(ztn, topt, xmax, zmax):
    zheight = np.zeros((zmax, xmax))
    for x in range(0, xmax):
        for z in range(0, zmax):
            zheight[z, x] = ztn[z] * \
                (1-topt[x]/ztn[zmax-1])+topt[x]

    return zheight


def pressure(pi):
    r"""Calculate the pressure from Exner function using PI."""
    p0 = 1000.
    cp = 1004.
    R = 287.

    return p0*np.power((pi/cp), (cp/R))


def temperature(theta, pi):
    r"""Calculate the temperature from Exner function using THETA and PI."""
    cp = 1004.

    return theta*(pi/cp)
