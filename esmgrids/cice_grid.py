
import numpy as np
import netCDF4 as nc

from .base_grid import BaseGrid

from os import environ
from datetime import datetime

from .util import *

def create_nc(filename):

    f = nc.Dataset(filename, 'w')

    f.timeGenerated = f"{datetime.now()}"
    f.created_by = f"{environ.get('USER')}"
    if is_git_repo():
        git_url, _ , git_hash = git_info()
        f.history = f"Created using commit {git_hash} of {git_url}"

    return f 

class CiceGrid(BaseGrid):

    def __init__(self, **kwargs):
        self.type = 'Arakawa B'
        self.full_name = 'CICE'

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

            dx_tn = f.variables['htn'][:] * 100.0
            dy_te = f.variables['hte'][:] * 100.0

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
                mask_t = f.variables['kmt'][:]

        return cls(x_t=x_t, y_t=y_t, x_u=x_u, y_u=y_u,
                   dx_t=dx_tn, dy_t=dy_te,
                   dx_tn=dx_tn, dy_te=dy_te,
                   area_t=area_t, area_u=area_u,
                   clat_t=clat_t, clon_t=clon_t, clat_u=clat_u, clon_u=clon_u,
                   mask_t=mask_t, description=description)

    def create_gridnc(self, grid_filename):
        self.grid_f = create_nc(grid_filename)
        return True
    
    def create_masknc(self, mask_filename):
        self.mask_f = create_nc(mask_filename)
        return True
    
    def write(self):
        """
        Write out CICE grid to netcdf.
        """

        f = self.grid_f

        # Create dimensions.
        f.createDimension('nx', self.num_lon_points)
        f.createDimension('ny', self.num_lat_points)

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

        # Convert from m to cm.
        htn[:] = self.dx_tn[:] * 100.
        hte[:] = self.dy_te[:] * 100.

        angle[:] = np.deg2rad(self.angle_u[:])
        angleT[:] = np.deg2rad(self.angle_t[:])

        f.close()

    def write_mask(self):

        """
        Write out CICE mask/kmt to netcdf.
        """

        f = self.mask_f
        f.createDimension('nx', self.num_lon_points)
        f.createDimension('ny', self.num_lat_points)
        mask = f.createVariable('kmt', 'f8', dimensions=('ny', 'nx'))
        # CICE uses 0 as masked, whereas internally we use 1 as masked.
        mask[:] = (1 - self.mask_t)
        f.close()


