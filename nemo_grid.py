#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from base_grid import BaseGrid

class NemoGrid(BaseGrid):

    def __init__(self, h_grid_def, v_grid_def=None, mask_file=None,
                    description='NEMO tripolar'):

        self.type = 'Arakawa C'

        with nc.Dataset(h_grid_def) as f:

            # Get t-points.
            x_t = f.variables['glamt'][:]
            y_t = f.variables['gphit'][:]

            x_u = f.variables['glamu'][:]
            y_u = f.variables['gphiu'][:]

            x_v = f.variables['glamv'][:]
            y_v = f.variables['gphiv'][:]

            # These variables hold the corners
            self.glamf = f.variables['glamf'][:]
            self.gphif = f.variables['gphif'][:]

            # These hold the edges.
            self.e1t = f.variables['e1t'][:]
            self.e2t = f.variables['e2t'][:]
            self.e1u = f.variables['e1u'][:]
            self.e2u = f.variables['e2u'][:]
            self.e1v = f.variables['e1v'][:]
            self.e2v = f.variables['e2v'][:]

        z = [0]
        if v_grid_def is not None:
            with nc.Dataset(v_grid_def) as f:
                z = f.variables['depth'][:]

        self.mask_file = mask_file

        super(NemoGrid, self).__init__(x_t, y_t, levels=z, x_u=x_u, y_u=y_u,
                x_v=x_v, y_v=y_v, description=description)

    def set_mask(self):
        """
        Mask is read from file if it exists, otherwise default to super class action.
        """

        if self.mask_file is None:
            super(NemoGrid, self).set_mask()
        else:
            with nc.Dataset(mask_file) as f:
                mask = np.zeros_like(f.variables['mask'], dtype=bool)
                mask[f.variables['mask'][:] == 0.0] = True
                self.mask_t = mask
                self.mask_u = mask
                self.mask_v = mask

    def calc_areas(self):

        self.area_t = self.e1t[:] * self.e2t[:]
        self.area_u = self.e1u[:] * self.e2u[:]
        self.area_v = self.e1v[:] * self.e2v[:]

    def make_corners(self):

        # These are the top righ-hand corner of t cells.
        glamf = self.glamf
        gphif = self.gphif

        # Extend south so that Southern most cells can have bottom corners.
        gphif_new = np.ndarray((gphif.shape[0] + 1, gphif.shape[1] + 1))
        gphif_new[1:, 1:] = gphif[:]
        gphif_new[0, 1:] = gphif[0, :] - abs(gphif[1, :] - gphif[0, :])

        glamf_new = np.ndarray((glamf.shape[0] + 1, glamf.shape[1] + 1))
        glamf_new[1:, 1:] = glamf[:]
        glamf_new[0, 1:] = glamf[0, :]

        # Repeat first longitude so that Western most cells have left corners.
        gphif_new[:, 0] = gphif_new[:, -1]
        glamf_new[:, 0] = glamf_new[:, -1]

        gphif = gphif_new
        glamf = glamf_new

        # Corners of t points. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clon[:] = np.NAN
        clon[:,:,0] = glamf[0:-1,0:-1]
        clon[:,:,1] = glamf[0:-1,1:]
        clon[:,:,2] = glamf[1:,1:]
        clon[:,:,3] = glamf[1:,0:-1]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = gphif[0:-1,0:-1]
        clat[:,:,1] = gphif[0:-1,1:]
        clat[:,:,2] = gphif[1:,1:]
        clat[:,:,3] = gphif[1:,0:-1]
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat

        # FIXME: do these.
        self.clon_u = clon
        self.clat_u = clat

        self.clon_v = clon
        self.clat_v = clat
