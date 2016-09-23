#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from base_grid import BaseGrid

class MomGrid(BaseGrid):

    def __init__(self, h_grid_def, v_grid_def, mask_file, description):
        """
        See src/mom5/ocean_core/ocean_grids.F90 and
        MOM4_guide.pdf for a description of the mosaic MOM5 grid.
        """

        with nc.Dataset(h_grid_def) as f:

            # Select points from double density horizontal grid. Only
            # need t-points.
            x_t = f.variables['x'][1::2,1::2]
            y_t = f.variables['y'][1::2,1::2]
            self.x_vt = f.variables['x'][:]
            self.y_vt = f.variables['y'][:]

        with nc.Dataset(v_grid_def) as f:
            # Only take cell centres.
            z = f.variables['zeta'][1::2]

        with nc.Dataset(mask_file) as f:
            mask = np.zeros_like(f.variables['mask'], dtype=bool)
            mask[f.variables['mask'][:] == 0.0] = True

        super(MomGrid, self).__init__(x_t, y_t, z, mask, description)

        self.num_lat_points = y_t.shape[0]
        self.num_lon_points = y_t.shape[1]

    def make_corners(self):

        # Uses double density grid to figure out corners.
        x = self.x_vt
        y = self.y_vt

        # Corners of t points. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clon[:] = np.NAN
        clon[:,:,0] = x[0:-1:2,0:-1:2]
        clon[:,:,1] = x[0:-1:2,2::2]
        clon[:,:,2] = x[2::2,2::2]
        clon[:,:,3] = x[2::2,0:-1:2]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = y[0:-1:2,0:-1:2]
        clat[:,:,1] = y[0:-1:2,2::2]
        clat[:,:,2] = y[2::2,2::2]
        clat[:,:,3] = y[2::2,0:-1:2]
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat

"""
Some old code.

class MOMGrid:

    def __init__(self, grid_filename, mask_filename, output_dir):
        """
        All the work gets done here. Other grids can then access any MOM field
        with <mom_object>.field
        """

        self.grid_filename = os.path.join(output_dir,
                                          os.path.basename(grid_filename))
        self.mask_filename = os.path.join(output_dir,
                                          os.path.basename(mask_filename))

        shutil.copyfile(mask_filename, self.mask_filename)
        shutil.copyfile(grid_filename, self.grid_filename)

        with nc.Dataset(mask_filename) as f:
            self.mask = np.copy(f.variables['mask'])

        with nc.Dataset(grid_filename) as f:
            self.dy = np.copy(f.variables['dy'])
            self.dx = np.copy(f.variables['dx'])
            self.angle_dx = np.copy(f.variables['angle_dx'])

            self.make_corners(f)
            self.calc_t_and_u_areas(f)

    def calc_t_and_u_areas(self, f):
        #
        #Calculate (add up) areas of T and U cells using the ocean areas. 
        #

        area = np.copy(f.variables['area'])
        self.area_t = np.zeros((area.shape[0]/2, area.shape[1]/2))
        self.area_u = np.zeros((area.shape[0]/2, area.shape[1]/2))

        # Add up areas, going clockwise from bottom left.
        self.area_t = area[0::2, 0::2] + area[1::2, 0::2] + \
                      area[1::2, 1::2] + area[0::2, 1::2]

        # These need to wrap around the globe. Copy ocn_area and add an extra
        # column at the end.
        area_ext = np.append(area[:], area[:, 0:1], axis=1)
        self.area_u = area_ext[0::2, 1::2] + area_ext[1::2, 1::2] + \
                      area_ext[1::2, 2::2] + area_ext[0::2, 2::2]


    def make_corners(self, f):
        # The standard mom grid includes t-cell corners be specifying the u, v
        # grid. Here we extract that and put it into the format expected by
        # the regridder and OASIS.

        x = np.copy(f.variables['x'])
        y = np.copy(f.variables['y'])

        self.clon = np.zeros((4, x.shape[0] / 2, x.shape[1] / 2))
        self.clon[:] = np.NAN
        self.clat = np.zeros((4, x.shape[0] / 2, x.shape[1] / 2))
        self.clat[:] = np.NAN

        # Corner lats. 0 is bottom left and then counter-clockwise. 
        # This is the OASIS convention. 
        self.clat[0,:,:] = y[0:-1:2,0:-1:2]
        self.clat[1,:,:] = y[0:-1:2,2::2]
        self.clat[2,:,:] = y[2::2,2::2]
        self.clat[3,:,:] = y[2::2,0:-1:2]

        # Corner lons.
        self.clon[0,:,:] = x[0:-1:2,0:-1:2]
        self.clon[1,:,:] = x[0:-1:2,2::2]
        self.clon[2,:,:] = x[2::2,2::2]
        self.clon[3,:,:] = x[2::2,0:-1:2]

        # Select points from double density grid. Southern most U points are
        # excluded. also the last (Eastern) U points, they are duplicates of
        # the first.
        self.x_t = x[1::2,1::2]
        self.y_t = y[1::2,1::2]
        self.x_u = x[2::2,0:-1:2]
        self.y_u = y[2::2,0:-1:2]
"""
