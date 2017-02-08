
from __future__ import print_function

import numpy as np
from base_grid import BaseGrid

# Various regirdding techniques (e.g. SCRIP) don't cope well with singularities.
SOUTHERN_EXTENT = -89.9995
NORTHERN_EXTENT = 89.9995

class RegularGrid(BaseGrid):

    def __init__(self, num_lons, num_lats, mask=None,
                 levels=[0], description=''):

        dx = 360.0 / num_lons
        dy = (abs(SOUTHERN_EXTENT) + abs(NORTHERN_EXTENT)) / (num_lats - 1)
        dx_half = dx / 2
        dy_half = dy / 2

        # Set lats and lons.
        lons = np.linspace(0, 360, num_lons, endpoint=False)
        # lat points exclude the poles.
        lats = np.linspace(SOUTHERN_EXTENT + dy_half, NORTHERN_EXTENT - dy_half, num_lats)

        super(RegularGrid, self).__init__(lons, lats, mask=mask,
                                          levels=levels, description=description)
