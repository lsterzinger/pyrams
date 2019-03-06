# from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def plot_heights(heights, lats, lons, dx, save, filename):
    '''Plots snowfall amount variable across domain with map background'''

    # Get number of x and y datapoints
    ny = int(lats.shape[0])
    nx = int(lats.shape[1])

    width = nx*dx
    height = ny*dx

    # Find center lat/lon
    centlat = lats[int(ny/2), int(nx/2)]
    centlon = lons[int(ny/2), int(nx/2)]
  
    # Create figure
    fig = plt.figure(figsize=(16, 9))

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
    height_plot = m.contour(x, y, heights)
    plt.clabel(height_plot, inline=True)
    if save is True:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


# def plot_snow(snow, lats, lons, dx):
    '''Plots snowfall amount variable across domain with map background'''

    snow = np.ma.masked_equal(snow, 0)

    # Get number of x and y datapoints
    ny = int(lats.shape[0])
    nx = int(lats.shape[1])

    width = nx*dx
    height = ny*dx

    # Find center lat/lon
    centlat = lats[int(ny/2), int(nx/2)]
    centlon = lons[int(ny/2), int(nx/2)]

    # Convert from kg/m^3 to inches
    #snow = snow*0.0393701

    nws_precip_colors = [
        "#04e9e7",  # 0.01 - 0.10 inches
        "#019ff4",  # 0.10 - 0.25 inches
        "#0300f4",  # 0.25 - 0.50 inches
        "#02fd02",  # 0.50 - 0.75 inches
        "#01c501",  # 0.75 - 1.00 inches
        "#008e00",  # 1.00 - 1.50 inchesbi  
        "#fdf802",  # 1.50 - 2.00 inchesbi  
        "#e5bc00",  # 2.00 - 2.50 inchesbi  
        "#fd9500",  # 2.50 - 3.00 inchesbi  
        "#fd0000",  # 3.00 - 4.00 inchesbi  
        "#d40000",  # 4.00 - 5.00 inchesbi  
        "#bc0000",  # 5.00 - 6.00 inches
        "#f800fd",  # 6.00 - 8.00 inches
        "#9854c6",  # 8.00 - 10.00 inches
        "#fdfdfd"   # 10.00+
    ]

    precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)

    # snow_levels = [0.01,
    #                0.10,
    #                0.25,
    #                0.5,
    #                0.75,
    #                1.0,
    #                1.5,
    #                2.0,
    #                2.5,
    #                3.0,
    #                4.0,
    #                5.0,
    #                6.0,
    #                8.0,
    #                10.0]

    snow_levels = range(0, 75, 5)
    # Create figure
    fig = plt.figure(figsize=(16, 16))

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
    snow_plot = m.contourf(x, y, snow, snow_levels, cmap=precip_colormap)
    cbar = plt.colorbar(snow_plot, label='Snow [inches]')
    return snow_plot


# def init_map(lats, lons, dx):
    
    ny = int(lats.shape[0])
    nx = int(lats.shape[1])

    width = nx*dx
    height = ny*dx

    # Find center lat/lon
    centlat = lats[int(ny/2), int(nx/2)]
    centlon = lons[int(ny/2), int(nx/2)]

    fig = plt.figure(figsize=(16, 16))

    m = Basemap(projection='stere', lon_0=centlon, lat_0=centlat,
                lat_ts=centlat, width=width, height=height)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    parallels = np.arange(0., 90, 10.)
    m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
    meridians = np.arange(180., 360, 10.)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)

    return m