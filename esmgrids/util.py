
import numpy as np
import pyproj
from shapely.geometry import shape

proj_str = '+proj=laea +lat_0={} +lon_0={} +ellps=sphere'

def calc_area_of_polygons(clons, clats):
    """
    Calculate the area of lat-lon polygons.

    We project sphere onto a flat surface using an equal area projection
    and then calculate the area of flat polygon.

    This is slow we should do some caching to avoid recomputing.
    """

    areas = np.zeros(clons.shape[1:])
    areas[:] = np.NAN

    for j in range(areas.shape[0]):
        for i in range(areas.shape[1]):

            lats = clats[:, j, i]
            lons = clons[:, j, i]

            lat_centre = lats[0] + abs(lats[2] - lats[1]) / 2
            lon_centre = lons[0] + abs(lons[1] - lons[0]) / 2

            pa = pyproj.Proj(proj_str.format(lat_centre, lon_centre))
            x, y = pa(lons, lats)

            cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
            areas[j, i] = shape(cop).area

    assert(np.sum(areas) is not np.NAN)
    assert(np.min(areas) > 0)

    return areas
