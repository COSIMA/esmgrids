#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
import exceptions
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import mpl_toolkits.basemap as basemap

EARTH_AREA = 510072000e6

class BaseGrid(object):

    def __init__(self, x_t, y_t, levels=[0], x_u=None, y_u=None, x_v=None, y_v=None,
                 mask=None, description=''):

        if len(x_t.shape) == 1:
            # We expect this to be a regular grid.
            assert np.allclose(np.diff(x_t),
                     np.array([x_t[1] - x_t[0]]*(len(x_t)-1)))
            assert np.allclose(np.diff(y_t),
                     np.array([y_t[1] - y_t[0]]*(len(y_t)-1)), atol=1e-4)
            # Turn into tiled
            self.x_t = np.tile(x_t, (y_t.shape[0], 1))
            self.y_t = np.tile(y_t, (x_t.shape[0], 1))
            self.y_t = self.y_t.transpose()

            self.dy = abs(self.y_t[1, 0] - self.y_t[0, 0])
            self.dx = abs(self.x_t[0, 1] - self.x_t[0, 0])
        else:
            # There is not constant dx, dy.
            self.dy = None
            self.dx = None

            self.x_t = x_t
            self.y_t = y_t

        self.num_lat_points = self.x_t.shape[0]
        self.num_lon_points = self.x_t.shape[1]
        self.num_levels = len(levels)

        self.x_u = x_u
        self.y_u = y_u
        self.x_v = x_v
        self.y_v = y_v
        self.z = levels
        self.mask = mask
        self.description = description

        self.set_mask()
        self.make_corners()
        self.calc_areas()

    def set_mask(self):

        if self.mask is None:
            # Default is all unmasked, up to user to mask.
            self.mask_t = np.zeros((self.num_levels,
                                 self.num_lat_points, self.num_lon_points),
                                 dtype='int')
        else:
            self.mask_t = self.mask

        # FIXME: is a correct generalisation?
        self.mask_u = self.mask_t
        self.mask_v = self.mask_t

    def calc_areas(self):
        """
        """

        self.area_t = self.calc_area_of_polygons(self.clon_t, self.clat_t)
        assert(abs(1 - np.sum(self.area_t) / EARTH_AREA) < 5e-4)

    def make_corners(self):

        x = self.x_t
        y = self.y_t

        dx_half = self.dx / 2.0
        dy_half = self.dy / 2.0

        # Set grid corners, we do these one corner at a time. Start at the 
        # bottom left and go anti-clockwise. This is the SCRIP convention.
        clon = np.empty((4, self.num_lat_points, self.num_lon_points))
        clon[:] = np.NAN
        clon[0,:,:] = x - dx_half
        clon[1,:,:] = x + dx_half
        clon[2,:,:] = x + dx_half
        clon[3,:,:] = x - dx_half
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((4, self.num_lat_points, self.num_lon_points))
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

        self.clon_t = clon
        self.clat_t = clat


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
        if len(self.mask.shape) == 2:
            imask[:] = np.invert(self.mask[:])
        else:
            imask[:] = np.invert(self.mask[0, :, :])

        corner_lat = f.createVariable('corner_lat', 'f8',
                                      ('lats', 'lons', 'grid_corners'))
        corner_lat.units = 'degrees'
        corner_lat[:] = clat[:]

        corner_lon = f.createVariable('corner_lon', 'f8',
                                      ('lats', 'lons', 'grid_corners'))
        corner_lon.units = 'degrees'
        corner_lon[:] = clon[:]

        f.close()


    def write_scrip(self, filename, mask=None, write_test_scrip=True, history=''):
        """
        Write out grid in SCRIP format.
        """

        self.make_corners()

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
            mask = self.mask

        if len(mask.shape) == 2:
            imask[:] = np.invert(mask[:]).flatten()
        else:
            imask[:] = np.invert(mask[0, :, :]).flatten()

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

        x, y = m(clons, clats)

        for j in range(x.shape[1]):
            for i in range(x.shape[2]):
                areas[j, i] = area_polygon(zip(x[:, j, i], y[:, j, i]))

        assert(np.sum(areas) is not np.NAN)
        assert(np.min(areas) > 0)

        return areas

