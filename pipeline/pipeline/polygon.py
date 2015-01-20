from shapely.geometry.polygon import Polygon
from shapely.geometry import Point, MultiPoint
import shapely.ops as ops
import numpy as np
from load import wrap_america
from geoloc2 import load_countries
from mpl_toolkits.basemap import Basemap


def get_polygon(polygon, wrap=True):
    """get_polygon
    gets the boundary polygon tu run eems on
    
    Parameters
    ----------
    polygon : int or str
        if str, file with polygon in eems format
        if set to `0`, get a  convex hull around the sample points
        if set to `1`, use the countries the points are in
        if set to `2`, get continent data
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
        polygon = polygon[..., ::-1]
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
        meta_data['POINTS'] = [Point(a[1]['LONG'], a[1]['LAT']) for a in
                               meta_data.iterrows()]
    return meta_data


def apply_map_projection(polygon, meta_data, projection):
    m = Basemap(llcrnrlon=-50.,llcrnrlat=-50.,urcrnrlon=340.,urcrnrlat=65.,\
                resolution='h',projection=projection,
                lat_0=40.,lon_0=-20.,lat_ts=20.)
    polygon = ops.transform(m, polygon)
    meta_data.POINTS = [ops.transform(m, p) for p in meta_data.POINTS]
    meta_data.LONG = [p.coords[0][0] for p in meta_data.POINTS]
    meta_data.LAT = [p.coords[0][1] for p in meta_data.POINTS]

    return polygon, meta_data


def get_eems_area(params, meta_data):
    """get_eems_area
    The main function to obtain the region for eems to run on
    
    Parameters
    ----------
    params : Parameters
        a parameters object
    
    Returns
    -------
    polygon : shapely.geometry.Polygon
        a polygon describing the plotting region
    meta_data : pd.DataFrame
        the potentially filtered meta_data object
    """

    if params.polygon is None and params.region is not None:
        poly1 = get_region_polygon(params.region, params.map,
                                   params.region_buffer, params.wrap)
    elif params.polygon is not None and params.region is None:
        poly1 = get_polygon(params.polygon, wrap=True)
    elif params.polygon is None and params.region is None:
        poly1 = None
    elif params.polygon is not None and params.region is not None:
        poly_region = get_region_polygon(params.region, params.map,
                                         params.region_buffer, params.wrap)
        poly_file = get_polygon(params.polygon, wrap=True)
        print type(poly_region), type(poly_file)
        poly1 = poly_region.intersection(poly_file)


    # hulls depend on trhe data
    if poly1 is not None:
        meta_data = filter_individuals_based_on_location(meta_data, poly1)

    hull = None
    if params.hull:
        hull = get_hull(meta_data)
        hull = hull.buffer(params.sample_buffer)
        if poly1 is None:
            poly1 = hull

    elif params.envelope:
        hull = get_envelope(meta_data)
        hull = hull.buffer(params.sample_buffer)
        if poly1 is None:
            poly1 = hull

    if hull is not None:
        poly1 = poly1.intersection(hull)

    if params.map_projection is not None:
        poly1, meta_data = apply_map_projection(poly1, meta_data, params.map_projection)

    return poly1, meta_data


def get_hull(meta_data):
    pts = MultiPoint(list(meta_data.POINTS))
    return pts.convex_hull


def get_envelope(meta_data):
    pts = MultiPoint(list(meta_data.POINTS))
    return pts.envelope


def get_region_polygon(region, map_file='', rbuffer=1, wrap=True):
    if region is None:
        return None

    countries = load_countries(map_file, wrap_americas=wrap)
    eems_region = countries[region]
    polygon = eems_region.get_boundary_polygon(min_area=0.9,
                                               buffer_lvl=rbuffer,
                                               return_type="polygon")
    return polygon
