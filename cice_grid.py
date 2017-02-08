
import numpy as np
import netCDF4 as nc

from base_grid import BaseGrid

class CiceGrid(BaseGrid):

    def __init__(self, **kwargs):
        super(CiceGrid, self).__init__(**kwargs)

    @classmethod
    def fromfile(cls, h_grid_def, mask_file=None,
                 description='CICE tripolar'):
        """
        Read in grid definition from file(s).
        """

        with nc.Dataset(h_grid_def) as f:
            x_t = np.rad2deg(f.variables['tlat'][:])
            y_t = np.rad2deg(f.variables['tlon'][:])

            x_u = np.rad2deg(f.variables['ulat'][:])
            y_u = np.rad2deg(f.variables['ulon'][:])

            dx_t = f.variables['htn'][:] * 100.0
            dy_t = f.variables['hte'][:] * 100.0

            area_t = f.variables['tarea'][:]
            area_u = f.variables['uarea'][:]

            angle_t = np.rad2deg(f.variables['angleT'][:])
            angle_u = np.rad2deg(f.variables['angle'][:])

            if 'clon_t' in f.variables:
                clon_t = f.variables['clon_t'][:]
                clat_t = f.variables['clat_t'][:]
                clon_u = f.variables['clon_u'][:]
                clat_u = f.variables['clat_u'][:]
            else:
                clon_t = clat_t = clon_u = clat_u = None

        if mask_file is not None:
            with nc.Dataset(mask_file) as f:
                mask_t = f.variables['mask'][:]

        return cls(x_t, y_t, x_u=x_u, y_u=y_u,
                   dx_t=dx_t, dy_t=dy_t,
                   area_t=area_t, area_u=area_u,
                   clat_t=clat_t, clon_t=clon_t, clat_u=clat_u, clon_u=clon_u,
                   mask_t=mask_t, description=description)


    def write(self, grid_filename, mask_filename):
        """
        Write out CICE grid to netcdf.
        """

        f = nc.Dataset(grid_filename, 'w')

        # Create dimensions.
        f.createDimension('nx', self.num_lon_points)
        f.createDimension('ny', self.num_lat_points)
        f.createDimension('nc', 4)

        # Make all CICE grid variables.
        ulat = f.createVariable('ulat', 'f8', dimensions=('ny', 'nx'))
        ulat.units = "radians"
        ulat.title = "Latitude of U points"
        ulon = f.createVariable('ulon', 'f8', dimensions=('ny', 'nx'))
        ulon.units = "radians"
        ulon.title = "Longitude of U points"
        tlat = f.createVariable('tlat', 'f8', dimensions=('ny', 'nx'))
        tlat.units = "radians"
        tlat.title = "Latitude of T points"
        tlon = f.createVariable('tlon', 'f8', dimensions=('ny', 'nx'))
        tlon.units = "radians"
        tlon.title = "Longitude of T points"

        if self.clon_t is not None:
            clon_t = f.createVariable('clon_t', 'f8',
                                      dimensions=('4', 'ny', 'nx'))
            clon_t.units = "radians"
            clon_t.title = "Longitude of T cell corners"
            clat_t = f.createVariable('clat_t', 'f8',
                                      dimensions=('4', 'ny', 'nx'))
            clat_t.units = "radians"
            clat_t.title = "Latitude of T cell corners"
            clon_u = f.createVariable('clon_u', 'f8',
                                      dimensions=('4', 'ny', 'nx'))
            clon_u.units = "radians"
            clon_u.title = "Longitude of U cell corners"
            clat_u = f.createVariable('clat_u', 'f8',
                                      dimensions=('4', 'ny', 'nx'))
            clat_u.units = "radians"
            clat_u.title = "Latitude of U cell corners"

        htn = f.createVariable('htn', 'f8', dimensions=('ny', 'nx'))
        htn.units = "cm"
        htn.title = "Width of T cells on North side."
        hte = f.createVariable('hte', 'f8', dimensions=('ny', 'nx'))
        hte.units = "cm"
        hte.title = "Width of T cells on East side."

        angle = f.createVariable('angle', 'f8', dimensions=('ny', 'nx'))
        angle.units = "radians"
        angle.title = "Rotation angle of U cells."
        angleT = f.createVariable('angleT', 'f8', dimensions=('ny', 'nx'))
        angleT.units = "radians"
        angleT.title = "Rotation angle of T cells."

        area_t = f.createVariable('tarea', 'f8', dimensions=('ny', 'nx'))
        area_t.units = "m^2"
        area_t.title = "Area of T cells."
        area_u = f.createVariable('uarea', 'f8', dimensions=('ny', 'nx'))
        area_u.units = "m^2"
        area_u.title = "Area of U cells."

        area_t[:] = self.area_t[:]
        area_u[:] = self.area_u[:]

        # Convert units: degrees -> radians.
        tlat[:] = np.deg2rad(self.y_t)
        tlon[:] = np.deg2rad(self.x_t)
        ulat[:] = np.deg2rad(self.y_u)
        ulon[:] = np.deg2rad(self.x_u)

        if self.clon_t is not None:
            clon_t[:] = np.deg2rad(self.clon_t)
            clat_t[:] = np.deg2rad(self.clat_t)
            clon_u[:] = np.deg2rad(self.clon_u)
            clat_u[:] = np.deg2rad(self.clat_u)

        # Convert from m to cm.
        htn[:] = self.dx_t * 100.
        hte[:] = self.dy_t * 100.

        angle[:] = np.deg2rad(self.angle_u[:])
        angleT[:] = np.deg2rad(self.angle_t[:])

        f.close()

        with nc.Dataset(mask_filename, 'w') as f:
            f.createDimension('nx', self.num_lon_points)
            f.createDimension('ny', self.num_lat_points)
            mask = f.createVariable('kmt', 'f8', dimensions=('ny', 'nx'))
            mask[:] = self.mask_t
