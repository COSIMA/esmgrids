import numpy as np
import netCDF4 as nc
from .base_grid import BaseGrid


class NemoGrid(BaseGrid):

    def __init__(self, h_grid_def, v_grid_def=None, mask_file=None, description="NEMO tripolar"):
        self.type = "Arakawa C"
        self.full_name = "NEMO"

        with nc.Dataset(h_grid_def) as f:

            # Get t-points.
            x_t = f.variables["glamt"][:]
            y_t = f.variables["gphit"][:]

            x_u = f.variables["glamu"][:]
            y_u = f.variables["gphiu"][:]

            x_v = f.variables["glamv"][:]
            y_v = f.variables["gphiv"][:]

            # These variables hold the corners
            x_f = f.variables["glamf"][:]
            y_f = f.variables["gphif"][:]

            # These hold the edges.
            e1t = f.variables["e1t"][:]
            e2t = f.variables["e2t"][:]
            e1u = f.variables["e1u"][:]
            e2u = f.variables["e2u"][:]
            e1v = f.variables["e1v"][:]
            e2v = f.variables["e2v"][:]

        z = [0]
        if v_grid_def is not None:
            with nc.Dataset(v_grid_def) as f:
                z = f.variables["depth"][:]

        with nc.Dataset(mask_file) as f:
            # NEMO default is 0 masked, 1 not masked.
            # Internal representation is True for masked, False
            # not masked.
            mask_t = 1 - f.variables["tmask"][:, :, :, :]
            mask_u = 1 - f.variables["umask"][:, :, :, :]
            mask_v = 1 - f.variables["vmask"][:, :, :, :]

        area_t = e1t[:] * e2t[:]
        area_u = e1u[:] * e2u[:]
        area_v = e1v[:] * e2v[:]

        clat_t, clon_t, clat_u, clon_u, clat_v, clon_v = make_corners(x_f, y_f, x_t, y_t, x_u, y_u, x_v, y_v)

        super(NemoGrid, self).__init__(
            x_t=x_t,
            y_t=y_t,
            x_u=x_u,
            y_u=y_u,
            x_v=x_v,
            y_v=y_v,
            area_t=area_t,
            area_u=area_u,
            area_v=area_v,
            clat_t=clat_t,
            clon_t=clon_t,
            clat_u=clat_u,
            clon_u=clon_u,
            clat_v=clat_v,
            clon_v=clon_v,
            mask_t=mask_t,
            mask_u=mask_u,
            mask_v=mask_v,
            levels=z,
            description=description,
        )


def make_corners(x_f, y_f, x_t, y_t, x_u, y_u, x_v, y_v):

    ##
    # Extend f points south so that south t cells can have bottom
    # corners. Also need to extend west to have left corners.
    ##
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

    y_u_new = np.ndarray((y_u.shape[0] + 1, y_u.shape[1] + 1))
    y_u_new[:-1, 1:] = y_u[:, :]
    y_u_new[-1, 1:] = y_u[-1, :]

    x_u_new = np.ndarray((x_u.shape[0] + 1, x_u.shape[1] + 1))
    x_u_new[:-1, 1:] = x_u[:, :]
    x_u_new[-1, 1:] = x_u[-1, :]

    # Repeat first longitude so that west t cells have left corners.
    y_u_new[:, 0] = y_u_new[:, -1]
    x_u_new[:, 0] = x_u_new[:, -1]

    y_u = y_u_new
    x_u = x_u_new

    # Corners of t cells are f points. Index 0 is bottom left and then
    # anti-clockwise.
    clon_t = np.empty((4, x_t.shape[0], x_t.shape[1]))
    clon_t[:] = np.NAN
    clon_t[0, :, :] = x_f[0:-1, 0:-1]
    clon_t[1, :, :] = x_f[0:-1, 1:]
    clon_t[2, :, :] = x_f[1:, 1:]
    clon_t[3, :, :] = x_f[1:, 0:-1]
    assert not np.isnan(np.sum(clon_t))

    clat_t = np.empty((4, x_t.shape[0], x_t.shape[1]))
    clat_t[:] = np.NAN
    clat_t[0, :, :] = y_f[0:-1, 0:-1]
    clat_t[1, :, :] = y_f[0:-1, 1:]
    clat_t[2, :, :] = y_f[1:, 1:]
    clat_t[3, :, :] = y_f[1:, 0:-1]
    assert not np.isnan(np.sum(clat_t))

    # The corners of u cells are v points.
    clon_u = np.empty((4, x_t.shape[0], x_t.shape[1]))
    clon_u[:] = np.NAN

    clon_u[0, :, :] = x_v[0:-1, 0:-1]
    clon_u[1, :, :] = x_v[0:-1, 1:]
    clon_u[2, :, :] = x_v[1:, 1:]
    clon_u[3, :, :] = x_v[1:, 0:-1]
    assert not np.isnan(np.sum(clon_u))

    clat_u = np.empty((4, x_t.shape[0], x_t.shape[1]))
    clat_u[:] = np.NAN
    clat_u[0, :, :] = y_v[0:-1, 0:-1]
    clat_u[1, :, :] = y_v[0:-1, 1:]
    clat_u[2, :, :] = y_v[1:, 1:]
    clat_u[3, :, :] = y_v[1:, 0:-1]
    assert not np.isnan(np.sum(clat_u))

    # The corners of v cells are u points.
    clon_v = np.empty((4, x_t.shape[0], x_t.shape[1]))
    clon_v[:] = np.NAN

    clon_v[0, :, :] = x_u[0:-1, 0:-1]
    clon_v[1, :, :] = x_u[0:-1, 1:]
    clon_v[2, :, :] = x_u[1:, 1:]
    clon_v[3, :, :] = x_u[1:, 0:-1]
    assert not np.isnan(np.sum(clon_v))

    clat_v = np.empty((4, x_t.shape[0], x_t.shape[1]))
    clat_v[:] = np.NAN
    clat_v[0, :, :] = y_u[0:-1, 0:-1]
    clat_v[1, :, :] = y_u[0:-1, 1:]
    clat_v[2, :, :] = y_u[1:, 1:]
    clat_v[3, :, :] = y_u[1:, 0:-1]
    assert not np.isnan(np.sum(clat_v))

    return clat_t, clon_t, clat_u, clon_u, clat_v, clon_v
