
import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid

class Jra55Grid(BaseGrid):

    def __init__(self, h_grid_def, description='JRA55 regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'JRA55'

        with nc.Dataset(h_grid_def) as f:
            lon_bnds = f.variables['lon_bnds'][:] 
            lat_bnds = f.variables['lat_bnds'][:]
            dx_t = lon_bnds[:, 1] - lon_bnds[:, 0]
            dy_t = lat_bnds[:, 1] - lat_bnds[:, 0]
            x_t = lon_bnds[:, 0] + dx_t[:] / 2
            y_t = lat_bnds[:, 0] + dy_t[:] / 2

        super(Jra55Grid, self).__init__(x_t=x_t, y_t=y_t, description=description)

    def fix_pole_holes(self):

        clat_copy = np.copy(self.clat_t)

        self.clat_t[2, -1,:] = 90.0
        self.clat_t[3, -1,:] = 90.0

        # Do South pole as well
        self.clat_t[0, 0,:] = -90.0
        self.clat_t[1, 0,:] = -90.0

