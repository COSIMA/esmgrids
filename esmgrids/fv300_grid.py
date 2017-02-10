#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from .regular_grid import RegularGrid

class FV300Grid(RegularGrid):

    def __init__(self, num_lons=128, num_lats=64, num_levels=1, mask_file=None, description=''):

        levels = range(num_levels)

        self.type = 'Arakawa A'
        self.full_name = 'FV300'

        if mask_file:
            with nc.Dataset(mask_file) as f:
                try:
                    mask_t = np.round(f.variables['WGOCN'][0, 0, :, :-1])
                except KeyError as e:
                    print("Error: var WGOCN not in {}.".format(mask_file),
                           file=sys.stderr)
                    raise e
        else:
            mask_t = np.ones((num_lats, num_lons))

        assert mask_t.shape[0] == num_lats
        assert mask_t.shape[1] == num_lons

        super(FV300Grid, self).__init__(num_lons, num_lats, mask_t=mask_t,
                                        levels=levels, description=description)
