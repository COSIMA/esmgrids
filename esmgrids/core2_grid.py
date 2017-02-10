
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
