
import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid


class Jra55RiverGrid(BaseGrid):

    def __init__(self, h_grid_def, description='JRA55 river regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'JRA55_river'

        with nc.Dataset(h_grid_def) as f:
            x_t = f.variables['longitude'][:]
            y_t = f.variables['latitude'][:]

        super(Jra55RiverGrid, self).__init__(x_t=x_t, y_t=y_t,
                                             description=description)

    def fix_pole_holes(self):

        self.clat_t[2, -1, :] = 90.0
        self.clat_t[3, -1, :] = 90.0

        # Do South pole as well
        self.clat_t[0, 0, :] = -90.0
        self.clat_t[1, 0, :] = -90.0
