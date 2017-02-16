
import numpy as np
import netCDF4 as nc
from .base_grid import BaseGrid

class Core2Grid(BaseGrid):

    def __init__(self, h_grid_def, description='CORE2 regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'CORE2'

        with nc.Dataset(h_grid_def) as f:
            x_t = f.variables['LON'][:]
            y_t = f.variables['LAT'][:]

        super(Core2Grid, self).__init__(x_t=x_t, y_t=y_t, description=description)


    def fix_pole_holes(self):
        """
        The CORE2 is not global, it has 'holes' at the poles. This will cause a
        problem when doing a conservative mapping to the MOM tripolar grid
        because the MOM grid includes the North Pole.

        The simplest way to deal with this is to extend the top corners of the
        northmost cells right to the pole. We do the same in the South.
        """

        clat_copy = np.copy(self.clat_t)

        self.clat_t[2, -1,:] = 90.0
        self.clat_t[3, -1,:] = 90.0

        # Do South pole as well
        self.clat_t[0, 0,:] = -90.0
        self.clat_t[1, 0,:] = -90.0

        assert ((self.clat_t[2, -1, :] - clat_copy[2, -1, :]) < 0.6).all()
        assert ((self.clat_t[3, -1, :] - clat_copy[3, -1, :]) < 0.6).all()
        assert ((self.clat_t[0, 0, :] - clat_copy[0, 0, :]) < 0.6).all()
        assert ((self.clat_t[1, 0, :] - clat_copy[1, 0, :]) < 0.6).all()
