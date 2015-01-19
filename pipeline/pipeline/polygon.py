from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import numpy as np
from load import wrap_america


def get_polygon(polygon, expand=0, wrap=True):
    """get_polygon
    gets the boundary polygon tu run eems on
    
    Parameters
    ----------
    polygon : int or str
        if str, file with polygon in eems format
        if set to `0`, get a  convex hull around the sample points
        if set to `1`, use the countries the points are in
        if set to `2`, get continent data
    expand : float
        how much the polygon should be expanded
    wrap : bool
        should coordinates in the americas be wrapped s.t they appear in the
        east?
    
    Returns
    -------
    polygon: shapely Polygon object
        the polygon to be used in the eems analysis
    """
    if polygon == 0:
        pass
    elif polygon == 1:
        raise NotImplementedError("country mode nyi")
    elif polygon == 2:
        raise NotImplementedError("continent mode nyi")
    else:
        polygon = np.loadtxt(polygon)
        polygon[:, 1] = wrap_america(polygon[:, 1])
        return Polygon(polygon)


def filter_individuals_based_on_location(meta_data, polygon):
    """filter_individuals
    only retains individuals that are inside the polygon
    
    Parameters
    ----------
    meta_data : pd.data_frame
        data frame with individuals and their location
    polygon : list of (int, int)
        a polygon, output from get_polygon
    
    Returns
    -------
    filtered_meta_data : pd.DataFrmae
        a data frame only with individuals used for analysis
    """
    create_points(meta_data)
    to_keep = [polygon.contains(pt) for pt in meta_data['POINTS']]
        
    return meta_data[to_keep]


def create_points(meta_data):
    """ adds a point object for meta data """
    if 'POINTS' not in meta_data:
        meta_data['POINTS'] = [Point(a[1]['LAT'], a[1]['LONG']) for a in
                               meta_data.iterrows()]
    return meta_data
