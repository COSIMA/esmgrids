
from __future__ import print_function

import numpy as np
import netCDF4 as nc
from base_grid import BaseGrid

class WoaGrid(BaseGrid):

    def __init__(self, grid_def):

        with nc.Dataset(grid_def) as f:
            x_t = f.variables['lon'][:]
            y_t = f.variables['lat'][:]
            depth = f.variables['depth'][:]
            mask = f.variables['t_an'][0, :, :, :].mask

        super(WoaGrid, self).__init__(x_t=x_t, y_t=y_t, levels=depth, mask=mask,
                                          description='WOA 1 degree grid')
