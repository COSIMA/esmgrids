
import numpy as np
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import mpl_toolkits.basemap as basemap

def calc_area_of_polygons(clons, clats):
    """
    Calculate the area of lat-lon polygons.

    We project sphere onto a flat surface using an equal area projection
    and then calculate the area of flat polygon.
    """

    def area_polygon(p):
        """
        Calculate the area of a polygon.

        Input is a polygon represented as a list of (x,y) vertex
        coordinates, implicitly wrapping around from the last vertex to the
        first.

        See http://stackoverflow.com/questions/451426/how-do-i-calculate-the-surface-area-of-a-2d-polygon
        """

        def segments(v):
            return zip(v, v[1:] + [v[0]])

        return 0.5 * abs(sum(x0*y1 - x1*y0
                             for ((x0, y0), (x1, y1)) in segments(p)))


    areas = np.zeros(clons.shape[1:])
    areas[:] = np.NAN

    m = basemap.Basemap(projection='laea', resolution='h',
                        llcrnrlon=0, llcrnrlat=-90.0,
                        urcrnrlon=360, urcrnrlat=90.0, lat_0=-90, lon_0=0)

    x, y = m(clons, clats)

    for j in range(x.shape[1]):
        for i in range(x.shape[2]):
            areas[j, i] = area_polygon(list(zip(x[:, j, i], y[:, j, i])))

    assert(np.sum(areas) is not np.NAN)
    assert(np.min(areas) > 0)

    return areas
