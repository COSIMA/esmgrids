#!/usr/bin/env python

from __future__ import print_function

import sys
import numpy as np
import netCDF4 as nc
from base_grid import BaseGrid

class T42Grid(BaseGrid):

    def __init__(self, num_lons, num_lats, num_levels, mask_file, description=''):

        # Set lats and lons.
        lons = np.linspace(0, 360, num_lons, endpoint=False)
        lats = np.linspace(-90, 90, num_lats)
        levels = range(num_levels)

        self.type = 'Spectral'

        with nc.Dataset(mask_file) as f:
            try:
                mask = np.round(f.variables['WGOCN'][:, :-1])
            except KeyError, e:
                print("Error: can't find ocean fraction var WGOCN in {}.".format(mask_file),
                       file=sys.stderr)
                raise e

        super(T42Grid, self).__init__(lons, lats, levels, mask, description)
