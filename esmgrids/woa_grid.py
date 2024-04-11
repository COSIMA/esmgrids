from __future__ import print_function

import numpy as np
import netCDF4 as nc
from .base_grid import BaseGrid


class WoaGrid(BaseGrid):

    def __init__(self, grid_def, calc_areas=True):

        with nc.Dataset(grid_def) as f:
            x_t = f.variables["lon"][:]
            y_t = f.variables["lat"][:]

            try:
                depth = f.variables["depth"][:]
            except (KeyError, ValueError) as e:
                depth = [0]

            try:
                mask = f.variables["t_an"][0, :, :, :].mask
            except (KeyError, ValueError) as e:
                try:
                    mask = f.variables["so"][0, :, :, :].mask
                except (KeyError, ValueError) as e:
                    mask = None

        super(WoaGrid, self).__init__(
            x_t=x_t, y_t=y_t, mask_t=mask, levels=depth, calc_areas=calc_areas, description="WOA 1 degree grid"
        )
