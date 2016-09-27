#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from base_grid import BaseGrid

class MomGrid(BaseGrid):

    def __init__(self, h_grid_def, v_grid_def, mask_file, description='MOM tripolar'):
        """
        See src/mom5/ocean_core/ocean_grids.F90 and
        MOM4_guide.pdf for a description of the mosaic MOM5 grid.
        """

        self.type = 'Arakawa B'

        with nc.Dataset(h_grid_def) as f:

            # Select points from double density horizontal grid.
            # t-cells.
            x_t = f.variables['x'][1::2,1::2]
            y_t = f.variables['y'][1::2,1::2]

            # u-cells.
            x_u = f.variables['x'][:-1:2,:-1:2]
            y_u = f.variables['y'][:-1:2,:-1:2]

            # All points
            self.x_dd = f.variables['x'][:]
            self.y_dd = f.variables['y'][:]

        with nc.Dataset(v_grid_def) as f:
            # Only take cell centres.
            z = f.variables['zeta'][1::2]

        with nc.Dataset(mask_file) as f:
            mask = np.zeros_like(f.variables['mask'], dtype=bool)
            mask[f.variables['mask'][:] == 0.0] = True

        super(MomGrid, self).__init__(x_t, y_t, z, x_u=x_u, y_u=y_u, mask=mask,
                                        description=description)

    def make_corners(self):

        print('Make corners on derived class called')

        # Uses double density grid to figure out corners.
        x = self.x_dd
        y = self.y_dd

        # Corners of t cells. Index 0 is bottom left and then
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

        # Corners of u cells. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_u.shape[0], self.x_u.shape[1], 4))
        clon[:] = np.NAN
        # The Southernmost row of u cells is half-size.
        # FIXME
        assert(False)
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

        self.clon_u = clon
        self.clat_u = clat

