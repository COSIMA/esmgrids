import numpy as np
import netCDF4 as nc
import exception

from base_grid import BaseGrid

class MomGrid(BaseGrid):
    """
    See src/mom5/ocean_core/ocean_grids.F90 and
    MOM4_guide.pdf for a description of the mosaic MOM5 grid.
    """

    def __init__(*args, **kwargs):
        super(MomGrid, self).__init__(args, kwargs)

    @classmethod
    def fromfile(self, h_grid_def, v_grid_def=None, mask_file=None,
                 description='MOM tripolar'):
        """
        Read in grid definition from file(s).
        """

        self.type = 'Arakawa B'

        with nc.Dataset(h_grid_def) as f:

            if 'x_T' in f.variables:
                # This is an old style horizontal grid definition.
                # t-cells.
                x_t = f.variables['x_T'][:]
                y_t = f.variables['y_T'][:]

                # u-cells.
                x_u = f.variables['x_C'][:]
                y_u = f.variables['y_C'][:]

                self.area_t = f.variables['area_T'][:]
                self.area_u = f.variables['area_C'][:]

                self.clon_t = f.variables['x_vert_T'][:]
                self.clat_t = f.variables['y_vert_T'][:]
                self.clon_u = f.variables['x_vert_C'][:]
                self.clat_u = f.variables['y_vert_C'][:]

            else:

                # Select points from double density horizontal grid.
                # t-cells.
                x_t = f.variables['x'][1::2,1::2]
                y_t = f.variables['y'][1::2,1::2]

                # u-cells.
                x_u = f.variables['x'][:-1:2,:-1:2]
                y_u = f.variables['y'][:-1:2,:-1:2]

                angle_t = f.variables['angle_dx'][1::2,1::2]
                angle_u = f.variables['angle_dx'][2::2,0:-1:2]

                area = f.variables['area'][:]
                self.area_t = np.zeros((area.shape[0]/2, area.shape[1]/2))
                self.area_u = np.zeros((area.shape[0]/2, area.shape[1]/2))

                # Add up areas, going clockwise from bottom left.
                self.area_t = area[0::2, 0::2] + area[1::2, 0::2] + \
                              area[1::2, 1::2] + area[0::2, 1::2]

                # These need to wrap around the globe. Copy ocn_area and
                # add an extra column at the end.
                area_ext = np.append(area[:], area[:, 0:1], axis=1)
                self.area_u = area_ext[0::2, 1::2] + area_ext[1::2, 1::2] + \
                              area_ext[1::2, 2::2] + area_ext[0::2, 2::2]

                x = f.variables['x'][:]
                y = f.variables['y'][:]
                clat_t, clon_t, clat_u, clon_u = make_corners(x, y)

        z = [0]
        if v_grid_def is not None:
            with nc.Dataset(v_grid_def) as f:
                # Only take cell centres.
                z = f.variables['zeta'][1::2]

        if mask_file is not None:
            with nc.Dataset(self.mask_file) as f:
                if 'wet' in f.variables:
                    self.mask_t = f.variables['wet'][:]
                    self.mask_u = f.variables['wet'][:]
                else:
                    mask = np.zeros_like(f.variables['mask'], dtype=bool)
                    mask[f.variables['mask'][:] >= 0.5] = True
                    self.mask_t = mask
                    self.mask_u = mask

        return cls(x_t, y_t, x_u=x_u, y_u=y_u, area_t=area_t, area_u=area_u,
                   clat_t=clat_t, clon_t=clon_t, clat_u=clat_u, clon_u=clon_u,
                   angle_t=angle_t, angle_u=angle_u,
                   mask_t=mask_t, mask_u=mask_u, levels=z,
                   description=description)

    def make_corners(self, x, y):

        # Corners of t cells. Index 0 is bottom left and then
        # anti-clockwise.
        clon_t = np.empty((4, self.x_t.shape[0], self.x_t.shape[1]))
        clon_t[:] = np.NAN
        clon_t[0,:,:] = x[0:-1:2,0:-1:2]
        clon_t[1,:,:] = x[0:-1:2,2::2]
        clon_t[2,:,:] = x[2::2,2::2]
        clon_t[3,:,:] = x[2::2,0:-1:2]
        assert(not np.isnan(np.sum(clon)))

        clat_t = np.empty((4, self.x_t.shape[0], self.x_t.shape[1]))
        clat_t[:] = np.NAN
        clat_t[0,:,:] = y[0:-1:2,0:-1:2]
        clat_t[1,:,:] = y[0:-1:2,2::2]
        clat_t[2,:,:] = y[2::2,2::2]
        clat_t[3,:,:] = y[2::2,0:-1:2]
        assert(not np.isnan(np.sum(clat_t)))

        # Corners of u cells. Index 0 is bottom left and then
        # anti-clockwise.

        # Need to be careful with the edges.
        # - Make the South most row of cells half size in the vertical.
        # - West needs to wrap around.
        # Do the easy bits first and then fix up below.

        clon_u = np.empty((4, self.x_u.shape[0], self.x_u.shape[1]))
        clon_u[:] = np.NAN
        clon_u[0,1:,1:] = x[1:-2:2,1:-2:2]
        clon_u[1,1:,1:] = x[1:-2:2,3::2]
        clon_u[2,1:,1:] = x[3::2,3::2]
        clon_u[3,1:,1:] = x[3::2,1:-2:2]

        # Fix up bottom row excluding left most column
        clon_u[0,0,1:] = x[0,1:-2:2]
        clon_u[1,0,1:] = x[0,3::2]
        clon_u[2,0,1:] = x[1,3::2]
        clon_u[3,0,1:] = x[1,1:-2:2]

        # Fix up leftmost column excluding bottom row
        clon_u[0,1:,0] = x[1:-2:2,-1]
        clon_u[1,1:,0] = x[1:-2:2,1]
        clon_u[2,1:,0] = x[3::2,1]
        clon_u[3,1:,0] = x[3::2,-1]

        # Fix up the bottom left corner point
        clon_u[0, 0, 0] = x[0,-1]
        clon_u[1, 0, 0] = x[0,1]
        clon_u[2, 0, 0] = x[1,1]
        clon_u[3, 0, 0] = x[1,-1]
        assert(not np.isnan(np.sum(clon_u)))

        clat_u = np.empty((4, self.x_t.shape[0], self.x_t.shape[1]))
        clat_u[:] = np.NAN
        clat_u[0,1:,1:] = y[1:-2:2,1:-2:2]
        clat_u[1,1:,1:] = y[1:-2:2,3::2]
        clat_u[2,1:,1:] = y[3::2,3::2]
        clat_u[3,1:,1:] = y[3::2,1:-2:2]

        # Fix up bottom row excluding left most column
        clat_u[0,0,1:] = y[0,1:-2:2]
        clat_u[1,0,1:] = y[0,3::2]
        clat_u[2,0,1:] = y[1,3::2]
        clat_u[3,0,1:] = y[1,1:-2:2]

        # Fix up leftmost column excluding bottom row
        clat_u[0,1:,0] = y[1:-2:2,-1]
        clat_u[1,1:,0] = y[1:-2:2,1]
        clat_u[2,1:,0] = y[3::2,1]
        clat_u[3,1:,0] = y[3::2,-1]

        # Fix up the bottom left corner point
        clat_u[0,0, 0] = y[0,-1]
        clat_u[1,0, 0] = y[0,1]
        clat_u[2,0, 0] = y[1,1]
        clat_u[3,0, 0] = y[1,-1]
        assert(not np.isnan(np.sum(clat_u)))

        return clat_t, clon_t, clat_u, clon_u

    def write(self, filename):
        """
        Write out a MOM grid to a netcdf file. The input can be any other
        grid type. For example this could be used to write out a MOM grid
        with the same characteristics as a CICE grid.
        """

        raise exception.NotImplementedError
