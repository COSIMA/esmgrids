
from __future__ import print_function

import sys
import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid

class DaitrenRunoffGrid(BaseGrid):

    def __init__(self, h_grid_def, description='Daitren runoff regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'Daitren_runoff'

        try:
            with nc.Dataset(h_grid_def) as f:
                x_t = f.variables['xc'][:]
                y_t = f.variables['yc'][:]
                clon_t = f.variables['xv'][:]
                clat_t = f.variables['yv'][:]
                mask_t = f.variables['mask'][:]
                area_t = f.variables['area'][:]
        except IOError:
            print('Error opening {}'.format(h_grid_def), file=sys.stderr)
            sys.exit(1)

        super(DaitrenRunoffGrid, self).__init__(x_t=x_t, y_t=y_t,
                                                clon_t=clon_t, clat_t=clat_t,
                                                mask_t=mask_t, area_t=area_t,
                                                description=description)
