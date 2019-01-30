import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid


class MomGrid(BaseGrid):
    """
    See src/mom5/ocean_core/ocean_grids.F90 and
    MOM4_guide.pdf for a description of the mosaic MOM5 grid.
    """

    def __init__(self, **kwargs):
        self.type = 'Arakawa B'
        self.full_name = 'MOM'

        super(MomGrid, self).__init__(**kwargs)

    @classmethod
    def fromfile(cls, h_grid_def, v_grid_def=None, mask_file=None,
                 description='MOM tripolar'):
        """
        Read in grid definition from file(s).
        """

        with nc.Dataset(h_grid_def) as f:

            if 'x_T' in f.variables:
                # This is an old style horizontal grid definition.
                # t-cells.
                x_t = f.variables['x_T'][:]
                y_t = f.variables['y_T'][:]

                # u-cells.
                x_u = f.variables['x_C'][:]
                y_u = f.variables['y_C'][:]

                area_t = f.variables['area_T'][:]
                area_u = f.variables['area_C'][:]

                clon_t = f.variables['x_vert_T'][:]
                clat_t = f.variables['y_vert_T'][:]
                clon_u = f.variables['x_vert_C'][:]
                clat_u = f.variables['y_vert_C'][:]

                dx_t = dy_t = dx_u = dy_u = None
                angle_t = angle_u = None

            else:

                x = f.variables['x'][:]
                y = f.variables['y'][:]

                # Select points from double density horizontal grid.
                # t-cell centre points
                x_t = x[1::2, 1::2]
                y_t = y[1::2, 1::2]

                # u-cell centre points
                x_u = x[2::2, 2::2]
                y_u = y[2::2, 2::2]

                dx = f.variables['dx'][:]
                dy = f.variables['dy'][:]

                # Through the centre of t cells
                dx_t = dx[1::2, ::2] + dx[1::2, 1::2]
                dy_t = dy[::2, 1::2] + dy[1::2, 1::2]

                # Along the N and E edges of t cells
                dx_tn = dx[2::2, ::2] + dx[2::2, 1::2]
                dy_te = dy[::2, 2::2] + dy[1::2, 2::2]

                # These need to wrap around the globe in order to do dx_u
                # add an extra column at the end.
                dx_ext = np.append(dx[:], dx[:, 0:1], axis=1)

                # The u-cells cross the tri-polar fold
                dy_ext = np.append(dy[:], dy[-1:, :], axis=0)

                # Through the centre of u cells
                dx_u = dx_ext[2::2, 1::2] + dx_ext[2::2, 2::2]
                dy_u = dy_ext[1::2, 2::2] + dy_ext[2::2, 2::2]

                # Along the N and E edges of u cells
                dx_un = dx_ext[3::2, 1::2] + dx_ext[3::2, 2::2]
                dy_ue = dy_ext[1::2, 3::2] + dy_ext[2::2, 3::2]

                angle_dx = f.variables['angle_dx'][:]
                # The angle of rotation is a corner cell centres and applies to
                # both t and u cells.
                angle_t = angle_dx[2::2, 2::2]
                angle_u = angle_dx[2::2, 2::2]

                area = f.variables['area'][:]

                # Add up areas, going clockwise from bottom left.
                area_t = area[0::2, 0::2] + area[1::2, 0::2] + \
                    area[1::2, 1::2] + area[0::2, 1::2]

                # These need to wrap around the globe. Copy ocn_area and
                # add an extra column at the end. Also u-cells cross the
                # tri-polar fold so add an extra row at the top.
                area_ext = np.append(area[:], area[:, 0:1], axis=1)
                area_ext = np.append(area_ext[:], area_ext[-1:, :], axis=0)

                area_u = area_ext[1::2, 1::2] + area_ext[2::2, 1::2] + \
                    area_ext[2::2, 2::2] + area_ext[1::2, 2::2]

                clat_t, clon_t, clat_u, clon_u, _, _ = make_corners(x, y)

        z = [0]
        if v_grid_def is not None:
            with nc.Dataset(v_grid_def) as f:
                # Only take cell centres.
                z = f.variables['zeta'][1::2]

        if mask_file is not None:
            with nc.Dataset(mask_file) as f:
                # MOM default is 0 masked, 1 not masked.
                # Internal representation is True for masked, False
                # not masked.
                if 'wet' in f.variables:
                    mask = np.ones_like(f.variables['wet'], dtype=bool)
                    mask[f.variables['wet'][:] >= 0.5] = False
                else:
                    mask = np.ones_like(f.variables['mask'], dtype=bool)
                    mask[f.variables['mask'][:] >= 0.5] = False
                mask_t = mask
                # FIXME: this is not correct.
                mask_u = mask
        else:
            mask_t = None
            mask_u = None

        return cls(x_t=x_t, y_t=y_t, x_u=x_u, y_u=y_u,
                   dx_t=dx_t, dy_t=dy_t, dx_u=dx_u, dy_u=dy_u,
                   dx_tn=dx_tn, dy_te=dy_te,
                   area_t=area_t, area_u=area_u,
                   clat_t=clat_t, clon_t=clon_t, clat_u=clat_u, clon_u=clon_u,
                   angle_t=angle_t, angle_u=angle_u,
                   mask_t=mask_t, mask_u=mask_u, levels=z,
                   description=description)

    def write(self, filename):
        """
        Write out a MOM grid to a netcdf file. The input can be any other
        grid type. For example this could be used to write out a MOM grid
        with the same characteristics as a CICE grid.
        """

        raise NotImplementedError


def make_corners(x, y):

    nrow = x.shape[0] // 2
    ncol = x.shape[1] // 2

    # Corners of t cells. Index 0 is bottom left and then
    # anti-clockwise.
    clon_t = np.empty((4, nrow, ncol))
    clon_t[:] = np.NAN
    clon_t[0, :, :] = x[0:-1:2, 0:-1:2]
    clon_t[1, :, :] = x[0:-1:2, 2::2]
    clon_t[2, :, :] = x[2::2, 2::2]
    clon_t[3, :, :] = x[2::2, 0:-1:2]
    assert(not np.isnan(np.sum(clon_t)))

    clat_t = np.empty((4, nrow, ncol))
    clat_t[:] = np.NAN
    clat_t[0, :, :] = y[0:-1:2, 0:-1:2]
    clat_t[1, :, :] = y[0:-1:2, 2::2]
    clat_t[2, :, :] = y[2::2, 2::2]
    clat_t[3, :, :] = y[2::2, 0:-1:2]
    assert(not np.isnan(np.sum(clat_t)))

    # Corners of u cells. Index 0 is bottom left and then
    # anti-clockwise.

    # Need to be careful with the edges.
    # - Make the South most row of cells half size in the vertical.
    # - West needs to wrap around.
    # Do the easy bits first and then fix up below.

    clon_u = np.empty((4, nrow, ncol))
    clon_u[:] = np.NAN
    clon_u[0, 1:, 1:] = x[1:-2:2, 1:-2:2]
    clon_u[1, 1:, 1:] = x[1:-2:2, 3::2]
    clon_u[2, 1:, 1:] = x[3::2, 3::2]
    clon_u[3, 1:, 1:] = x[3::2, 1:-2:2]

    # Fix up bottom row excluding left most column
    clon_u[0, 0, 1:] = x[0, 1:-2:2]
    clon_u[1, 0, 1:] = x[0, 3::2]
    clon_u[2, 0, 1:] = x[1, 3::2]
    clon_u[3, 0, 1:] = x[1, 1:-2:2]

    # Fix up leftmost column excluding bottom row
    clon_u[0, 1:, 0] = x[1:-2:2, -1]
    clon_u[1, 1:, 0] = x[1:-2:2, 1]
    clon_u[2, 1:, 0] = x[3::2, 1]
    clon_u[3, 1:, 0] = x[3::2, -1]

    # Fix up the bottom left corner point
    clon_u[0, 0, 0] = x[0, -1]
    clon_u[1, 0, 0] = x[0, 1]
    clon_u[2, 0, 0] = x[1, 1]
    clon_u[3, 0, 0] = x[1, -1]
    assert(not np.isnan(np.sum(clon_u)))

    clat_u = np.empty((4, nrow, ncol))
    clat_u[:] = np.NAN
    clat_u[0, 1:, 1:] = y[1:-2:2, 1:-2:2]
    clat_u[1, 1:, 1:] = y[1:-2:2, 3::2]
    clat_u[2, 1:, 1:] = y[3::2, 3::2]
    clat_u[3, 1:, 1:] = y[3::2, 1:-2:2]

    # Fix up bottom row excluding left most column
    clat_u[0, 0, 1:] = y[0, 1:-2:2]
    clat_u[1, 0, 1:] = y[0, 3::2]
    clat_u[2, 0, 1:] = y[1, 3::2]
    clat_u[3, 0, 1:] = y[1, 1:-2:2]

    # Fix up leftmost column excluding bottom row
    clat_u[0, 1:, 0] = y[1:-2:2, -1]
    clat_u[1, 1:, 0] = y[1:-2:2, 1]
    clat_u[2, 1:, 0] = y[3::2, 1]
    clat_u[3, 1:, 0] = y[3::2, -1]

    # Fix up the bottom left corner point
    clat_u[0, 0,  0] = y[0, -1]
    clat_u[1, 0,  0] = y[0, 1]
    clat_u[2, 0,  0] = y[1, 1]
    clat_u[3, 0,  0] = y[1, -1]
    assert(not np.isnan(np.sum(clat_u)))

    return clat_t, clon_t, clat_u, clon_u, None, None
