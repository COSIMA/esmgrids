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

    def __init__(self, x_t, y_t, x_u=None, y_u=None, x_v=None, y_v=None,
                 dx_t=None, dy_t=None, dx_u=None, dy_u=None, dx_v=None, dy_v=None,
                 area_t=None, area_u=None, area_v=None,
                 clat_t=None, clon_t=None, clat_u=None, clon_u=None, 
                 clat_v=None, clon_v=None,
                 mask_t=None, mask_u=None, mask_v=None,
                 levels=[0], description=''):

        self.x_t = x_t; self.y_t = y_t; self.x_u = x_u; self.y_u = y_u
        self.x_v = x_v; self.y_v = y_v
        self.dx_t = dx_t; self.dy_t = dy_t
        self.dx_u = dx_u; self.dy_u = dy_u
        self.dx_v = dx_v; self.dy_v = dy_v
        self.area_t = area_t; self.area_u = area_u; self.area_v = area_v
        self.clat_t = clat_t; self.clon_t = clon_t
        self.clat_u = clat_u; self.clon_u = clon_u
        self.clat_v = clat_v; self.clon_v = clon_v
        self.mask_t = mask_t; self.mask_u = mask_u; self.mask_v = mask_v

        self.z = levels
        self.mask = mask
        self.description = description

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

            if not dx_t:
                self.dx_t = abs(self.x_t[0, 1] - self.x_t[0, 0])
            if not dy_t
                self.dy_t = abs(self.y_t[1, 0] - self.y_t[0, 0])

        if not self.mask_t:
            # Default is all unmasked, up to user to mask.
            self.mask_t = np.zeros((self.num_levels,
                                 self.num_lat_points, self.num_lon_points),
                                 dtype='int')
        if not self.mask_u:
            # FIXME
            self.mask_u = self.mask_t
        if not self.mask_v:
            # FIXME
            self.mask_v = self.mask_t

        self.num_lat_points = self.x_t.shape[0]
        self.num_lon_points = self.x_t.shape[1]
        self.num_levels = len(levels)

        if not self.clon_t or not self.clat_t:
            self.make_corners()

        if not self.area_t:
            self.calc_areas()

    def calc_areas(self):

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
                areas[j, i] = area_polygon(zip(x[:, j, i], y[:, j, i]))

        assert(np.sum(areas) is not np.NAN)
        assert(np.min(areas) > 0)

        return areas

