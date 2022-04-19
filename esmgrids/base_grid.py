
import numpy as np
import netCDF4 as nc

from .util import calc_area_of_polygons


class BaseGrid(object):

    def __init__(self, **kwargs):

        x_t = kwargs.get('x_t', None)
        y_t = kwargs.get('y_t', None)
        self.x_u = kwargs.get('x_u', None)
        self.y_u = kwargs.get('y_u', None)
        self.x_v = kwargs.get('x_v', None)
        self.y_v = kwargs.get('y_v', None)

        dx_t = kwargs.get('dx_t', None)
        dy_t = kwargs.get('dy_t', None)

        self.dx_tn = kwargs.get('dx_tn', None)
        self.dy_te = kwargs.get('dy_te', None)

        self.dx_u = kwargs.get('dx_u', None)
        self.dy_u = kwargs.get('dy_u', None)
        self.dx_v = kwargs.get('dx_v', None)
        self.dy_v = kwargs.get('dy_v', None)

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
        self.calc_areas = kwargs.get('calc_areas', True)

        if len(x_t.shape) == 1:
            # Tile x_t
            self.x_t = np.tile(x_t, (y_t.shape[0], 1))
            self.y_t = np.tile(y_t, (x_t.shape[0], 1))
            self.y_t = self.y_t.transpose()
        else:
            self.x_t = x_t
            self.y_t = y_t

        if dx_t is not None and len(dx_t.shape) == 1:
            # Tile dx_t
            self.dx_t = np.tile(x_t, (dy_t.shape[0], 1))
            self.dy_t = np.tile(y_t, (dx_t.shape[0], 1))
            self.dy_t = self.dy_t.transpose()
        elif dx_t is None:
            # Repeat the first row/column using the values closest.
            self.dx_t = np.ndarray((self.x_t.shape[0], self.x_t.shape[1]-1))
            self.dy_t = np.ndarray((self.x_t.shape[0]-1, self.x_t.shape[1]))
            self.dx_t[:, :] = self.x_t[:, 1:] - self.x_t[:, :-1]
            self.dy_t[:, :] = self.y_t[1:, :] - self.y_t[:-1, :]
        else:
            self.dx_t = dx_t
            self.dy_t = dy_t

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

        # The base class may implement this if it has a hole at the north
        # pole to be fixed.
        self.fix_pole_holes()

        # The base class may implement this if the corner latitudes overrun the
        # poles
        self.fix_pole_overruns()

        if self.area_t is None and self.calc_areas:
            self.area_t = calc_area_of_polygons(self.clon_t, self.clat_t)

    def fix_pole_holes(self):
        pass


    def fix_pole_overruns(self):
        pass


    @classmethod
    def fromgrid(cls, grid):
        """
        Read in definition from another grid-like object.
        """

        return cls(**grid.__dict__)

    def write_scrip(self, filename, mask=None,
                    write_test_scrip=True, history=''):
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
        corner_lat[:, 0] = clat[0, :, :].flatten()
        corner_lat[:, 1] = clat[1, :, :].flatten()
        corner_lat[:, 2] = clat[2, :, :].flatten()
        corner_lat[:, 3] = clat[3, :, :].flatten()

        corner_lon = f.createVariable('grid_corner_lon', 'f8',
                                      ('grid_size', 'grid_corners'))
        corner_lon.units = 'degrees'
        corner_lon[:, 0] = clon[0, :, :].flatten()
        corner_lon[:, 1] = clon[1, :, :].flatten()
        corner_lon[:, 2] = clon[2, :, :].flatten()
        corner_lon[:, 3] = clon[3, :, :].flatten()

        f.title = self.description
        f.history = history
        f.close()

        if write_test_scrip:
            self.write_test_scrip(filename + '_test')

    def write_test_scrip(self, filename):
        """
        Write out SCRIP grid contents in a format which is easier to read/test.
        """

        f = nc.Dataset(filename, 'w')

        x = self.x_t
        y = self.y_t
        clat = self.clat_t
        clon = self.clon_t

        f.createDimension('lats', self.num_lat_points)
        f.createDimension('lons', self.num_lon_points)
        f.createDimension('grid_corners', 4)
        f.createDimension('grid_rank', 2)

        center_lat = f.createVariable('center_lat', 'f8', ('lats', 'lons'))
        center_lat.units = 'degrees'
        center_lat[:] = y[:]

        center_lon = f.createVariable('center_lon', 'f8', ('lats', 'lons'))
        center_lon.units = 'degrees'
        center_lon[:] = x[:]

        imask = f.createVariable('mask', 'i4', ('lats', 'lons'))
        imask.units = 'unitless'
        # Invert the mask. SCRIP uses zero for points that do not
        # participate.
        if len(self.mask_t.shape) == 2:
            imask[:] = 1 - self.mask_t[:]
        else:
            imask[:] = 1 - self.mask_t[0, :, :]

        corner_lat = f.createVariable('corner_lat', 'f8',
                                      ('lats', 'lons', 'grid_corners'))
        corner_lat.units = 'degrees'
        corner_lat[:, :, 0] = clat[0, :, :].flatten()
        corner_lat[:, :, 1] = clat[1, :, :].flatten()
        corner_lat[:, :, 2] = clat[2, :, :].flatten()
        corner_lat[:, :, 3] = clat[3, :, :].flatten()

        corner_lon = f.createVariable('corner_lon', 'f8',
                                      ('lats', 'lons', 'grid_corners'))
        corner_lon.units = 'degrees'
        corner_lon[:, :, 0] = clon[0, :, :].flatten()
        corner_lon[:, :, 1] = clon[1, :, :].flatten()
        corner_lon[:, :, 2] = clon[2, :, :].flatten()
        corner_lon[:, :, 3] = clon[3, :, :].flatten()

        f.close()



def make_corners(x, y, dx, dy):

    dx_half = dx / 2.0
    dy_half = dy / 2.0
    nrow = x.shape[0]
    ncol = x.shape[1]

    # Set grid corners, we do these one corner at a time. Start at the
    # bottom left and go anti-clockwise. This is the SCRIP convention.
    clon = np.empty((4, nrow, ncol))
    clon[:] = np.NAN

    clon[0, :, 1:] = x[:, 1:] - dx_half[:, :]
    clon[0, :, 0] = x[:, 0] - dx_half[:, 0]
    clon[3, :, 1:] = x[:, 1:] - dx_half[:, :]
    clon[3, :, 0] = x[:, 0] - dx_half[:, 0]

    clon[1, :, :-1] = x[:, :-1] + dx_half[:, :]
    clon[1, :, -1] = x[:, -1] + dx_half[:, -1]
    clon[2, :, :-1] = x[:, :-1] + dx_half[:, :]
    clon[2, :, -1] = x[:, -1] + dx_half[:, -1]

    assert(not np.isnan(np.sum(clon)))

    clat = np.empty((4, nrow, ncol))
    clat[:] = np.NAN

    clat[0, 1:, :] = y[1:, :] - dy_half[:, :]
    clat[0, 0, :] = y[0, :] - dy_half[0, :]
    clat[1, 1:, :] = y[1:, :] - dy_half[:, :]
    clat[1, 0, :] = y[0, :] - dy_half[0, :]

    clat[2, :-1, :] = y[:-1, :] + dy_half[:, :]
    clat[2, -1, :] = y[-1, :] + dy_half[-1, :]
    clat[3, :-1, :] = y[:-1, :] + dy_half[:, :]
    clat[3, -1, :] = y[-1, :] + dy_half[-1, :]

    assert(not np.isnan(np.sum(clat)))

    # The bottom latitude band should always be Southern extent.
#    assert(np.all(clat[0, 0, :] == np.min(y) - dy_half[0, :]))
#    assert(np.all(clat[1, 0, :] == np.min(y) - dy_half[0, :]))

    # The top latitude band should always be Northern extent.
#    assert(np.all(clat[2, -1, :] == np.max(y) + dy_half[-1, :]))
#    assert(np.all(clat[3, -1, :] == np.max(y) + dy_half[-1, :]))

    return clat, clon, None, None, None, None
