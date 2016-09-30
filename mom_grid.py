#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from base_grid import BaseGrid

class MomGrid(BaseGrid):

    def __init__(self, h_grid_def, v_grid_def=None, mask_file=None, description='MOM tripolar'):
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

            self.area_t = f.variables['area'][1::2, 1::2]
            self.area_u = f.variables['area'][:-1:2, :-1:2]

        z = [0]
        if v_grid_def is not None:
            with nc.Dataset(v_grid_def) as f:
                # Only take cell centres.
                z = f.variables['zeta'][1::2]

        self.mask_file = mask_file

        super(MomGrid, self).__init__(x_t, y_t, z, x_u=x_u, y_u=y_u,
                                        description=description)

    def set_mask(self):
        """
        Mask is read from file if it exists, otherwise default to super class action.
        """

        if self.mask_file is None:
            super(MOMGrid, self).set_mask()
        else:
            with nc.Dataset(self.mask_file) as f:
                mask = np.zeros_like(f.variables['mask'], dtype=bool)
                mask[f.variables['mask'][:] == 0.0] = True
                self.mask_t = mask
                self.mask_u = mask


    def calc_areas(self):
        """
        Areas are read from file above.
        """
        pass

    def make_corners(self):

        # Uses double density grid to figure out corners.
        x = self.x_dd
        y = self.y_dd

        # Corners of t cells. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clon[:] = np.NAN
        clon[:,:,0] = x[:-1:2,:-1:2]
        clon[:,:,1] = x[:-1:2,2::2]
        clon[:,:,2] = x[2::2,2::2]
        clon[:,:,3] = x[2::2,:-1:2]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = y[:-1:2,:-1:2]
        clat[:,:,1] = y[:-1:2,2::2]
        clat[:,:,2] = y[2::2,2::2]
        clat[:,:,3] = y[2::2,:-1:2]
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat

        # Corners of u cells. Index 0 is bottom left and then
        # anti-clockwise.

        # Need to be careful with the edges.
        # - Make the South most row of cells half size in the vertical.
        # - West needs to wrap around.
        # Do the easy bits first and then fix up below.

        clon = np.empty((self.x_u.shape[0], self.x_u.shape[1], 4))
        clon[:] = np.NAN
        clon[1:,1:,0] = x[1:-2:2,1:-2:2]
        clon[1:,1:,1] = x[1:-2:2,3::2]
        clon[1:,1:,2] = x[3::2,3::2]
        clon[1:,1:,3] = x[3::2,1:-2:2]

        # Fix up bottom row excluding left most column
        clon[0,1:,0] = x[0,1:-2:2]
        clon[0,1:,1] = x[0,3::2]
        clon[0,1:,2] = x[1,3::2]
        clon[0,1:,3] = x[1,1:-2:2]

        # Fix up leftmost column excluding bottom row
        clon[1:,0,0] = x[1:-2:2,-1]
        clon[1:,0,1] = x[1:-2:2,1]
        clon[1:,0,2] = x[3::2,1]
        clon[1:,0,3] = x[3::2,-1]

        # Fix up the bottom left corner point
        clon[0, 0, 0] = x[0,-1]
        clon[0, 0, 1] = x[0,1]
        clon[0, 0, 2] = x[1,1]
        clon[0, 0, 3] = x[1,-1]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[1:,1:,0] = y[1:-2:2,1:-2:2]
        clat[1:,1:,1] = y[1:-2:2,3::2]
        clat[1:,1:,2] = y[3::2,3::2]
        clat[1:,1:,3] = y[3::2,1:-2:2]

        # Fix up bottom row excluding left most column
        clat[0,1:,0] = y[0,1:-2:2]
        clat[0,1:,1] = y[0,3::2]
        clat[0,1:,2] = y[1,3::2]
        clat[0,1:,3] = y[1,1:-2:2]

        # Fix up leftmost column excluding bottom row
        clat[1:,0,0] = y[1:-2:2,-1]
        clat[1:,0,1] = y[1:-2:2,1]
        clat[1:,0,2] = y[3::2,1]
        clat[1:,0,3] = y[3::2,-1]

        # Fix up the bottom left corner point
        clat[0, 0, 0] = y[0,-1]
        clat[0, 0, 1] = y[0,1]
        clat[0, 0, 2] = y[1,1]
        clat[0, 0, 3] = y[1,-1]
        assert(not np.isnan(np.sum(clat)))

        self.clon_u = clon
        self.clat_u = clat
