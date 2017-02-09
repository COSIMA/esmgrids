
import numpy as np
import netCDF4 as nc
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import mpl_toolkits.basemap as basemap

EARTH_AREA = 510072000e6

class BaseGrid(object):

    def __init__(self, **kwargs):

        x_t = kwargs.get('x_t', None); y_t = kwargs.get('y_t', None)
        self.x_u = kwargs.get('x_u', None); self.y_u = kwargs.get('y_u', None)
        self.x_v = kwargs.get('x_v', None); self.y_v = kwargs.get('y_v', None)

        self.dx_t = kwargs.get('dx_t', None); self.dy_t = kwargs.get('dy_t', None)
        self.dx_u = kwargs.get('dx_u', None); self.dy_u = kwargs.get('dy_u', None)
        self.dx_v = kwargs.get('dx_v', None); self.dy_v = kwargs.get('dy_v', None)

        self.area_t = kwargs.get('area_t', None)
        self.area_u = kwargs.get('area_u', None)
        self.area_v = kwargs.get('area_v', None)

        self.clat_t = kwargs.get('clat_t', None)
        self.clon_t = kwargs.get('clon_t', None)
        self.clat_u = kwargs.get('clat_u', None)
        self.clon_u = kwargs.get('clon_u', None)
        self.clat_v = kwargs.get('clat_v', None)
        self.clon_v = kwargs.get('clon_v', None)

        self.angle_t = kwargs.get('angle_t', None)
        self.angle_u = kwargs.get('angle_u', None)
        self.angle_v = kwargs.get('angle_v', None)

        self.mask_t = kwargs.get('mask_t', None)
        self.mask_u = kwargs.get('mask_u', None)
        self.mask_v = kwargs.get('mask_v', None)

        self.z = kwargs.get('levels', [0])
        self.description = kwargs.get('description', '')

        if len(x_t.shape) == 1:
            # FIXME this shouldn't be here, it should be in
            # regular_grid.py
            # We expect this to be a regular grid.
            assert np.allclose(np.diff(x_t),
                     np.array([x_t[1] - x_t[0]]*(len(x_t)-1)))
            assert np.allclose(np.diff(y_t),
                     np.array([y_t[1] - y_t[0]]*(len(y_t)-1)), atol=1e-4)
            # Turn into tiled
            self.x_t = np.tile(x_t, (y_t.shape[0], 1))
            self.y_t = np.tile(y_t, (x_t.shape[0], 1))
            self.y_t = self.y_t.transpose()

            if not self.dx_t:
                self.dx_t = abs(self.x_t[0, 1] - self.x_t[0, 0])
            if not self.dy_t:
                self.dy_t = abs(self.y_t[1, 0] - self.y_t[0, 0])
        else:
            self.x_t = x_t
            self.y_t = x_t

        self.num_lat_points = self.x_t.shape[0]
        self.num_lon_points = self.x_t.shape[1]
        self.num_levels = len(self.z)

        if self.mask_t is None:
            # Default is all unmasked, up to user to mask.
            self.mask_t = np.zeros((self.num_levels,
                                 self.num_lat_points, self.num_lon_points),
                                 dtype='int')
        if self.mask_u is None:
            # FIXME
            self.mask_u = self.mask_t
        if self.mask_v is None:
            # FIXME
            self.mask_v = self.mask_t

        if self.clon_t is None or self.clat_t is None:
            self.clat_t, self.clon_t, _, _, _, _ = \
                make_corners(self.x_t, self.y_t, self.dx_t, self.dy_t)

        if self.area_t is None:
            self.calc_areas()

    @classmethod
    def fromgrid(cls, grid):
        """
        Read in definition from another grid-like object.
        """

        return cls(**grid.__dict__)


    def calc_areas(self):

        self.area_t = self.calc_area_of_polygons(self.clon_t, self.clat_t)
        assert(abs(1 - np.sum(self.area_t) / EARTH_AREA) < 5e-4)


    def write_scrip(self, filename, mask=None, write_test_scrip=True, history=''):
        """
        Write out grid in SCRIP format.
        """

        f = nc.Dataset(filename, 'w')

        x = self.x_t
        y = self.y_t

        clat = self.clat_t
        clon = self.clon_t
        num_points = self.num_lat_points * self.num_lon_points

        f.createDimension('grid_size', num_points)
        f.createDimension('grid_corners', 4)
        f.createDimension('grid_rank', 2)

        grid_dims = f.createVariable('grid_dims', 'i4', ('grid_rank'))
        # SCRIP likes lon, lat
        grid_dims[:] = [self.num_lon_points, self.num_lat_points]

        center_lat = f.createVariable('grid_center_lat', 'f8', ('grid_size'))
        center_lat.units = 'degrees'
        center_lat[:] = y[:].flatten()

        center_lon = f.createVariable('grid_center_lon', 'f8', ('grid_size'))
        center_lon.units = 'degrees'
        center_lon[:] = x[:].flatten()

        imask = f.createVariable('grid_imask', 'i4', ('grid_size'))
        imask.units = 'unitless'

        # Invert the mask. SCRIP uses zero for points that do not
        # participate.
        if mask is not None:
            mask = mask
        else:
            mask = self.mask_t

        if len(mask.shape) == 2:
            imask[:] = 1 - mask[:].flatten()
        else:
            imask[:] = 1 - mask[0, :, :].flatten()

        corner_lat = f.createVariable('grid_corner_lat', 'f8',
                                      ('grid_size', 'grid_corners'))
        corner_lat.units = 'degrees'
        corner_lat[:] = clat[:].flatten()

        corner_lon = f.createVariable('grid_corner_lon', 'f8',
                                      ('grid_size', 'grid_corners'))
        corner_lon.units = 'degrees'
        corner_lon[:] = clon[:].flatten()

        f.title = self.description
        f.history = history
        f.close()

        if write_test_scrip:
            self.write_test_scrip(filename + '_test')


    def calc_area_of_polygons(self, clons, clats):
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

        # FIXME, don't do this. 
        clats[np.where(clats == 90.0)] = 89.9995

        x, y = m(clons, clats)

        for j in range(x.shape[1]):
            for i in range(x.shape[2]):
                areas[j, i] = area_polygon(list(zip(x[:, j, i], y[:, j, i])))

        assert(np.sum(areas) is not np.NAN)
        assert(np.min(areas) > 0)

        return areas

def make_corners(x, y, dx, dy):
    """
    FIXME: this should be in regular_grid.py
    """

    dx_half = dx / 2.0
    dy_half = dy / 2.0
    nrow = x.shape[0]
    ncol = x.shape[1]

    # Set grid corners, we do these one corner at a time. Start at the 
    # bottom left and go anti-clockwise. This is the SCRIP convention.
    clon = np.empty((4, nrow, ncol))
    clon[:] = np.NAN
    clon[0,:,:] = x - dx_half
    clon[1,:,:] = x + dx_half
    clon[2,:,:] = x + dx_half
    clon[3,:,:] = x - dx_half
    assert(not np.isnan(np.sum(clon)))

    clat = np.empty((4, nrow, ncol))
    clat[:] = np.NAN
    clat[0,:,:] = y - dy_half
    clat[1,:,:] = y - dy_half
    clat[2,:,:] = y + dy_half
    clat[3,:,:] = y + dy_half
    assert(not np.isnan(np.sum(clat)))

    # The bottom latitude band should always be Southern extent.
    assert(np.all(clat[0, 0, :] == np.min(y) - dy_half))
    assert(np.all(clat[1, 0, :] == np.min(y) - dy_half))

    # The top latitude band should always be Northern extent.
    assert(np.all(clat[2, -1, :] == np.max(y) + dy_half))
    assert(np.all(clat[3, -1, :] == np.max(y) + dy_half))

    return clat, clon, None, None, None, None

