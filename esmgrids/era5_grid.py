
import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid


class Era5Grid(BaseGrid):

    def __init__(self, h_grid_def, description='ERA5 regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'ERA5'

        with nc.Dataset(h_grid_def) as f:
            x_t = f.variables['longitude'][:]
            # ERA5 puts latitudes the wrong way around
            y_t = f.variables['latitude'][:]

        super(Era5Grid, self).__init__(x_t=x_t, y_t=y_t, calc_areas=False,
                                        description=description)

    def fix_pole_overruns(self):
        """
        ERA5 has grid cell centres at latitudes of -90, 90.
        """

        self.clat_t[np.where(self.clat_t > 90)] = 90
        self.clat_t[np.where(self.clat_t < -90)] = -90
