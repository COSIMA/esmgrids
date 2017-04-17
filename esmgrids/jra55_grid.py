
import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid

class Jra55Grid(BaseGrid):

    def __init__(self, h_grid_def, description='JRA55 regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'JRA55'

        with nc.Dataset(h_grid_def) as f:
            x_t = f.variables['lon'][:]
            y_t = f.variables['lat'][1:-1]

        super(Jra55Grid, self).__init__(x_t=x_t, y_t=y_t, description=description)

    def fix_pole_holes(self):

        clat_copy = np.copy(self.clat_t)

        self.clat_t[2, -1,:] = 90.0
        self.clat_t[3, -1,:] = 90.0

        # Do South pole as well
        self.clat_t[0, 0,:] = -90.0
        self.clat_t[1, 0,:] = -90.0

