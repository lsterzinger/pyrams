""" Contains tools to aide in plotting of RAMS data """
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shapely.geometry as sgeom
from copy import copy

def _find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.

    Taken from https://gist.github.com/ajdawson/dd536f786741e987ae4e
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """
    Draw ticks on the bottom x-axis of a Lambert Conformal projection.

    Taken from https://gist.github.com/ajdawson/dd536f786741e987ae4e

    Arguments
    ---------
    ax : matplotlib.axis
        Axis on which to add ticks

    ticks : list
        x-ticks to add
    """
    
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])
    

def lambert_yticks(ax, ticks):
    """
    Draw ticks on the left y-axis of a Lambert Conformal projection.
    
    Taken from https://gist.github.com/ajdawson/dd536f786741e987ae4e

    Arguments
    ---------
    ax : matplotlib.axis
        Axis on which to add ticks

    ticks : list
        y-ticks to add
    """    
    
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """
    Get the tick locations and labels for an axis of a Lambert Conformal projection.
    
    Taken from https://gist.github.com/ajdawson/dd536f786741e987ae4e
    """
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = _find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels