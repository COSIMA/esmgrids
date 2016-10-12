#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from regular_grid import RegularGrid

class FV300Grid(RegularGrid):

    def __init__(self, num_lons, num_lats, num_levels, mask_file, description=''):

        levels = range(num_levels)

        self.type = 'Arakawa A'

        with nc.Dataset(mask_file) as f:
            try:
                mask = np.round(f.variables['WGOCN'][0, 0, :, :])
            except KeyError, e:
                print("Error: can't find ocean fraction var WGOCN in {}.".format(mask_file),
                       file=sys.stderr)
                raise e

        super(FV300Grid, self).__init__(num_lons, num_lats, levels, mask=mask,
                                        description=description)
