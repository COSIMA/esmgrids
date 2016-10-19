#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from base_grid import BaseGrid

class Mom1Grid(BaseGrid):

    def __init__(self, h_grid_def, v_grid_def=None, mask_file=None,
                    description='MOM tripolar'):
        """
        MOM 1 degree grid.
        """

        self.type = 'Arakawa B'

        with nc.Dataset(h_grid_def) as f:

            # Select points from double density horizontal grid.
            # t-cells.
            x_t = f.variables['x_T'][:]
            y_t = f.variables['y_T'][:]

            # u-cells.
            x_u = f.variables['x_C'][:]
            y_u = f.variables['y_C'][:]

            self.area_t = f.variables['area_T'][:]
            self.area_u = f.variables['area_C'][:]

            self.mask_t = 1 - f.variables['wet'][:]
            self.mask_u = 1 - f.variables['wet'][:]

            self.clon_t = f.variables['x_vert_T'][:]
            self.clat_t = f.variables['y_vert_T'][:]
            self.clon_u = f.variables['x_vert_C'][:]
            self.clat_u = f.variables['y_vert_C'][:]

        super(Mom1Grid, self).__init__(x_t, y_t, levels=[0], x_u=x_u, y_u=y_u,
                                        description=description)

    def set_mask(self):
        """
        Mask is read from file above.
        """
        pass

    def calc_areas(self):
        """
        Areas are read from file above.
        """
        pass

    def make_corners(self):
        """
        Corners are read from file above.
        """
