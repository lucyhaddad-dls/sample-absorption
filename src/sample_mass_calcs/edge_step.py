from .formulas import _edges
from numpy.lib.stride_tricks import sliding_window_view
from scipy.interpolate import interp1d
from xraylib import EdgeEnergy
import numpy as np

import matplotlib.pyplot as plt

def find_edges(y:np.ndarray)->tuple[np.ndarray, np.ndarray]:
    """
    Find coordinates for all potential edges detected
    in y array.

    Arguments:
        y (np.ndarray): Array of f2/absorption coefficient values.

    Returns:
        tuple (tuple): tuple containing:
            where_edge (np.ndarray): Array of coordinates for edges detected.
            where_pre (np.ndarray): Array of coordinates for pre-edge region.

    Example:
        ```
        from xraydb import f2_chantler
        x = np.arange(600, 9000, 1)
        y = f2_chantler("Ni", x)
        edges, pre_edges = find_edges(y)
        ```
    """

    dy = np.gradient(y)
  

    filt = lambda y: np.array([0 if i < 1 else 1 for i in y])
    
    dny = np.concat(([0], np.diff(filt(dy))))

   
    where_edge = np.where(dny < 0)[0]
    where_pre = np.where(dny > 0)[0]

    return where_edge, where_pre


def get_bounds(x:np.ndarray,
               y:np.ndarray,
               edge_energy:float|int)\
                ->tuple[int, tuple[int, int], int]:
    """
    Find closest edge in absorption coefficient given\
    a provided edge energy.

    Arguments:
        x (np.ndarray): Energy array.
        y (np.ndarray): Photoabsorption values array.
        edge_energy (float): Energy at which absorption edge should be.

    Returns:
        tuple (tuple): tuple containing:
            edge (int): Closest coordinate in y where the edge is.
            bounds (tuple[int, int]): Pre- to post-edge bounds around the edge.
            pre (int): Pre-edge upper bound.

    """
    edges, pres = find_edges(y)

    edge_points = np.concat(([0], edges, [len(x)-1]))
    edge_energies = np.array([x[p] for p in edge_points])

    winview = sliding_window_view(edge_energies, window_shape=2)
    midpoints = [np.where(x <= e)[0][-1] for e in winview.mean(axis=1)]

    midpoint_arr = np.concat(([0], midpoints, [len(x)-1]))

    # find closest pre-edge point to absorption edge:

    diffs_pre = [(c, np.abs(x[c]-edge_energy)) for c in pres]
    diffs_pre = sorted(diffs_pre, key=lambda c:c[1])
    pre = diffs_pre[0][0]
    # find closest edge to indexed pre-edge:
    diffs_edge = [(c, x[c]-x[pre]) for c in edges]
    # take only positive differences
    diffs_edge = [d for d in diffs_edge if d[-1] > 0]
    diffs_edge = sorted(diffs_edge, key = lambda c:c[-1])
    edge = diffs_edge[0][0]

    # find two midpoints closest to edge:
    edge_no = np.where(edges >= edge)[0][0]
    bounds = [midpoint_arr[edge_no+1], midpoint_arr[edge_no+2]]
    if bounds[0] < 0: bounds[0] = 0
    if bounds[-1] > len(x) - 1: bounds[-1] = len(x) - 1
    if np.abs(bounds[0] - pre) <= 3:
        # change to be a bit more smart? look at gradient.
        bounds[0] -= 5

    return edge, bounds, pre

def fit_bkg(x:np.ndarray,
            y:np.ndarray,
            bounds:tuple[int, int],
            pre_idx:int,
             edge_idx:float, include:bool)->tuple[float, float, float]:
    """
    Fit a quadratic background to the section of absorption values\
    containing an absorption edge and calculate step height.

    Arguments:
        x (np.ndarray): Energy axis.
        y (np.ndarray): Absorption coefficient array.
        bounds (tuple[int, int]): Full span for pre- and post-edge region.
        pre_idx (int): Where to end the pre-edge region from.
        edge_idx (int): Index of where energy = absorbing atom's edge energy.
        include (bool): Whether to calculate the jump for the given atom.

    Returns:
        tuple (tuple): tuple containing:
            absjump (float): Calculated jump height in f2.
            absmax (float): Value of f2 at maximum about the range.
            abspre (float): Value of f2 at the pre-edge end.
    """
    lower = bounds[0]; upper = bounds[-1]
    pre_y = y[lower:pre_idx]; post_y = y[edge_idx:upper]
    pre_x = x[lower:pre_idx]
    max_idx = np.where(post_y == np.max(post_y))[0][0] + edge_idx
    absmax = y[max_idx]; abspre = y[pre_idx]

    if include == False:
        absjump = 0
    else:
        lo_idx = np.where(pre_y == np.min(pre_y))[0][0] + lower
        if lo_idx != lower:
            pre_y = y[lower:lo_idx]; pre_x = x[lower:lo_idx]
        
        bkg = interp1d(pre_x, pre_y, kind="quadratic", fill_value="extrapolate")(x)
        absjump = y[max_idx] - bkg[max_idx]

    return absjump, absmax, abspre

def check_multi_edges(x:np.ndarray,
                      y:np.ndarray, 
                      abs_idx:int)->bool:
    """
    Check if any other elements in a compound will have an absorption
    edge within +/- 50 eV of the absorbing atom.

    Arguments:
        x (np.ndarray): Energy axis.
        y (np.ndarray): Photoabsorption values.
        abs_idx (int): Index where the energy = edge energy.

    Returns:
        include (bool): Whether to include the element in further edge-step
                         calculations.
    """
    where_edge, where_pre = find_edges(y)
    abs_energy = x[abs_idx]
    diffs = [np.abs(x[e]-abs_energy) for e in where_edge]
    if np.min(diffs) <= 50:
        include = True
    else:
        include = False
    return include

def get_edge_vals(x:np.ndarray,
                  y:np.ndarray,
                  x_target:float|int,
                  include:bool)->tuple[float, float, float]:
    
    edge, bounds, pre = get_bounds(x, y, x_target)

    jump_height, ymax, ymin = fit_bkg(x, y, bounds,
                                      pre, edge,
                                      include)

    return jump_height, ymax, ymin
    
