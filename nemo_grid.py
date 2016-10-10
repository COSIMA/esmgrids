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
            self.x_f = f.variables['glamf'][:]
            self.y_f = f.variables['gphif'][:]

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
            with nc.Dataset(self.mask_file) as f:
                self.mask_t = ~(f.variables['tmask'][:, :, :, :])
                self.mask_u = ~(f.variables['umask'][:, :, :, :])
                self.mask_v = ~(f.variables['vmask'][:, :, :, :])

    def calc_areas(self):

        self.area_t = self.e1t[:] * self.e2t[:]
        self.area_u = self.e1u[:] * self.e2u[:]
        self.area_v = self.e1v[:] * self.e2v[:]

    def make_corners(self):

        ##
        # Extend f points south so that south t cells can have bottom
        # corners. Also need to extend west to have left corners.
        ##
        x_f = self.x_f
        y_f = self.y_f

        y_f_new = np.ndarray((y_f.shape[0] + 1, y_f.shape[1] + 1))
        y_f_new[1:, 1:] = y_f[:]
        y_f_new[0, 1:] = y_f[0, :]

        x_f_new = np.ndarray((x_f.shape[0] + 1, x_f.shape[1] + 1))
        x_f_new[1:, 1:] = x_f[:]
        x_f_new[0, 1:] = x_f[0, :]

        # Repeat first longitude so that west t cells have left corners.
        y_f_new[:, 0] = y_f_new[:, -1]
        x_f_new[:, 0] = x_f_new[:, -1]

        y_f = y_f_new
        x_f = x_f_new

        ##
        # Extend v points south so that south u cells can have bottom
        # corners. Also need to extend east to have right corners.
        ##
        x_v = self.x_v
        y_v = self.y_v

        y_v_new = np.ndarray((y_v.shape[0] + 1, y_v.shape[1] + 1))
        y_v_new[1:, :-1] = y_v[:]
        y_v_new[0, :-1] = y_v[0, :]

        x_v_new = np.ndarray((x_v.shape[0] + 1, x_v.shape[1] + 1))
        x_v_new[1:, :-1] = x_v[:]
        x_v_new[0, :-1] = x_v[0, :]

        # Repeat last longitude so that east cells have right corners.
        y_v_new[:, -1] = y_v_new[:, 0]
        x_v_new[:, -1] = x_v_new[:, 0]

        y_v = y_v_new
        x_v = x_v_new

        ##
        # Extend u points north so that north v cells can have top
        # corners. Also need to extend west to have left corners.
        ##
        x_u = self.x_u
        y_u = self.y_u

        y_u_new = np.ndarray((y_u.shape[0] + 1, y_u.shape[1] + 1))
        y_u_new[:-1, 1:] = y_u[:, :]
        y_u_new[-1, 1:] = y_u[-1, :]

        x_u_new = np.ndarray((x_u.shape[0] + 1, x_u.shape[1] + 1))
        y_u_new[:-1, 1:] = y_u[:, :]
        y_u_new[-1, 1:] = y_u[-1, :]

        # Repeat first longitude so that west t cells have left corners.
        y_u_new[:, 0] = y_u_new[:, -1]
        x_u_new[:, 0] = x_u_new[:, -1]

        y_u = y_u_new
        x_u = x_u_new

        # Corners of t cells are f points. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clon[:] = np.NAN
        clon[:,:,0] = x_f[0:-1,0:-1]
        clon[:,:,1] = x_f[0:-1,1:]
        clon[:,:,2] = x_f[1:,1:]
        clon[:,:,3] = x_f[1:,0:-1]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = y_f[0:-1,0:-1]
        clat[:,:,1] = y_f[0:-1,1:]
        clat[:,:,2] = y_f[1:,1:]
        clat[:,:,3] = y_f[1:,0:-1]
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat

        # The corners of u cells are v points.
        clon = np.empty((self.x_u.shape[0], self.x_u.shape[1], 4))
        clon[:] = np.NAN

        clon[:,:,0] = x_v[0:-1,0:-1]
        clon[:,:,1] = x_v[0:-1,1:]
        clon[:,:,2] = x_v[1:,1:]
        clon[:,:,3] = x_v[1:,0:-1]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_u.shape[0], self.x_u.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = y_v[0:-1,0:-1]
        clat[:,:,1] = y_v[0:-1,1:]
        clat[:,:,2] = y_v[1:,1:]
        clat[:,:,3] = y_v[1:,0:-1]
        assert(not np.isnan(np.sum(clat)))

        self.clon_u = clon
        self.clat_u = clat

        # The corners of u cells are v points.
        clon = np.empty((self.x_v.shape[0], self.x_v.shape[1], 4))
        clon[:] = np.NAN

        clon[:,:,0] = x_u[0:-1,0:-1]
        clon[:,:,1] = x_u[0:-1,1:]
        clon[:,:,2] = x_u[1:,1:]
        clon[:,:,3] = x_u[1:,0:-1]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_v.shape[0], self.x_v.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = y_u[0:-1,0:-1]
        clat[:,:,1] = y_u[0:-1,1:]
        clat[:,:,2] = y_u[1:,1:]
        clat[:,:,3] = y_u[1:,0:-1]
        assert(not np.isnan(np.sum(clat)))

        self.clon_v = clon
        self.clat_v = clat
