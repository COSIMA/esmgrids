import numpy as np
import netCDF4 as nc

from esmgrids.base_grid import BaseGrid
from esmgrids.mom_grid import MomGrid
from esmgrids.util import *


class CiceGrid(BaseGrid):

    def __init__(self, **kwargs):
        self.type = "Arakawa B"
        self.full_name = "CICE"

        super(CiceGrid, self).__init__(**kwargs)

    @classmethod
    def fromfile(cls, h_grid_def, mask_file=None, description="CICE tripolar"):
        """
        Read in grid definition from file(s).
        """

        with nc.Dataset(h_grid_def) as f:
            x_t = np.rad2deg(f.variables["tlat"][:])
            y_t = np.rad2deg(f.variables["tlon"][:])

            x_u = np.rad2deg(f.variables["ulat"][:])
            y_u = np.rad2deg(f.variables["ulon"][:])

            dx_tn = f.variables["htn"][:] * 100.0
            dy_te = f.variables["hte"][:] * 100.0

            area_t = f.variables["tarea"][:]
            area_u = f.variables["uarea"][:]

            angle_t = np.rad2deg(f.variables["angleT"][:])
            angle_u = np.rad2deg(f.variables["angle"][:])

            if "clon_t" in f.variables:
                clon_t = f.variables["clon_t"][:]
                clat_t = f.variables["clat_t"][:]
                clon_u = f.variables["clon_u"][:]
                clat_u = f.variables["clat_u"][:]
            else:
                clon_t = clat_t = clon_u = clat_u = None

        if mask_file is not None:
            with nc.Dataset(mask_file) as f:
                mask_t = f.variables["kmt"][:]

        return cls(
            x_t=x_t,
            y_t=y_t,
            x_u=x_u,
            y_u=y_u,
            dx_t=dx_tn,
            dy_t=dy_te,
            dx_tn=dx_tn,
            dy_te=dy_te,
            area_t=area_t,
            area_u=area_u,
            clat_t=clat_t,
            clon_t=clon_t,
            clat_u=clat_u,
            clon_u=clon_u,
            mask_t=mask_t,
            description=description,
        )

    def create_gridnc(self, grid_filename):
        self.grid_f = create_nc(grid_filename)
        return True

    def create_masknc(self, mask_filename):
        self.mask_f = create_nc(mask_filename)
        return True

    def create_2d_grid_var(self, name):
        # set chunksizes based on OM2 config
        # To-do: load these from a configuration file?
        if self.num_lon_points == 360:  # 1deg
            chunksizes = (150, 180)
        elif self.num_lon_points == 1440:  # 0.25deg
            chunksizes = (540, 720)
        elif self.num_lon_points == 3600:  # 0.01deg
            chunksizes = (270, 360)
        else:
            chunksizes = None

        return self.grid_f.createVariable(
            name,
            "f8",
            dimensions=("ny", "nx"),
            compression="zlib",
            complevel=1,
            chunksizes=chunksizes,
        )

    def write(self):
        """
        Write out CICE grid to netcdf.
        """

        f = self.grid_f

        # Create dimensions.
        f.createDimension("nx", self.num_lon_points)
        # nx is the grid_longitude but doesn't have a value other than its index
        f.createDimension("ny", self.num_lat_points)
        # ny is the grid_latitude but doesn't have a value other than its index

        # Make all CICE grid variables.
        # names are based on https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
        f.Conventions = "CF-1.6"

        ulat = self.create_2d_grid_var("ulat")
        ulat.units = "radians"
        ulat.long_name = "Latitude of U points"
        ulat.standard_name = "latitude"
        ulon = self.create_2d_grid_var("ulon")
        ulon.units = "radians"
        ulon.long_name = "Longitude of U points"
        ulon.standard_name = "longitude"
        tlat = self.create_2d_grid_var("tlat")
        tlat.units = "radians"
        tlat.long_name = "Latitude of T points"
        tlat.standard_name = "latitude"
        tlon = self.create_2d_grid_var("tlon")
        tlon.units = "radians"
        tlon.long_name = "Longitude of T points"
        tlon.standard_name = "longitude"

        htn = self.create_2d_grid_var("htn")
        htn.units = "cm"
        htn.long_name = "Width of T cells on North side."
        htn.coordinates = "ulat tlon"
        htn.grid_mapping = "crs"
        hte = self.create_2d_grid_var("hte")
        hte.units = "cm"
        hte.long_name = "Width of T cells on East side."
        hte.coordinates = "tlat ulon"
        hte.grid_mapping = "crs"

        angle = self.create_2d_grid_var("angle")
        angle.units = "radians"
        angle.long_name = "Rotation angle of U cells."
        angle.standard_name = "angle_of_rotation_from_east_to_x"
        angle.coordinates = "ulat ulon"
        angle.grid_mapping = "crs"
        angleT = self.create_2d_grid_var("angleT")
        angleT.units = "radians"
        angleT.long_name = "Rotation angle of T cells."
        angleT.standard_name = "angle_of_rotation_from_east_to_x"
        angleT.coordinates = "tlat tlon"
        angleT.grid_mapping = "crs"

        area_t = self.create_2d_grid_var("tarea")
        area_t.units = "m^2"
        area_t.long_name = "Area of T cells."
        area_t.standard_name = "cell_area"
        area_t.coordinates = "tlat tlon"
        area_t.grid_mapping = "crs"
        area_u = self.create_2d_grid_var("uarea")
        area_u.units = "m^2"
        area_u.long_name = "Area of U cells."
        area_u.standard_name = "cell_area"
        area_u.coordinates = "ulat ulon"
        area_u.grid_mapping = "crs"

        area_t[:] = self.area_t[:]
        area_u[:] = self.area_u[:]

        # Convert units: degrees -> radians.
        tlat[:] = np.deg2rad(self.y_t)
        tlon[:] = np.deg2rad(self.x_t)
        ulat[:] = np.deg2rad(self.y_u)
        ulon[:] = np.deg2rad(self.x_u)

        # Convert from m to cm.
        htn[:] = self.dx_tn[:] * 100.0
        hte[:] = self.dy_te[:] * 100.0

        angle[:] = np.deg2rad(self.angle_u[:])
        angleT[:] = np.deg2rad(self.angle_t[:])

        f.close()

    def write_mask(self):
        """
        Write out CICE mask/kmt to netcdf.
        """

        f = self.mask_f
        f.createDimension("nx", self.num_lon_points)
        f.createDimension("ny", self.num_lat_points)
        mask = f.createVariable("kmt", "f8", dimensions=("ny", "nx"), compression="zlib", complevel=1)

        mask.grid_mapping = "crs"
        mask.standard_name = "sea_binary_mask"

        # CICE uses 0 as masked, whereas internally we use 1 as masked.
        mask[:] = 1 - self.mask_t
        f.close()


def cice_from_mom(ocean_hgrid, ocean_mask, grid_file="grid.nc", mask_file="kmt.nc"):

    mom = MomGrid.fromfile(ocean_hgrid, mask_file=ocean_mask)

    cice = CiceGrid.fromgrid(mom)

    cice.create_gridnc(grid_file)

    # Add versioning information
    cice.grid_f.inputfile = f"{ocean_hgrid}"
    cice.grid_f.inputfile_md5 = md5sum(ocean_hgrid)
    cice.grid_f.history_command = f"python make_CICE_grid.py {ocean_hgrid} {ocean_mask}"

    # Add the typical crs (i.e. WGS84/EPSG4326 , but in radians).
    crs = cice.grid_f.createVariable("crs", "S1")
    crs.grid_mapping_name = "tripolar_latitude_longitude"
    crs.crs_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["radians",1,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

    cice.write()

    cice.create_masknc(mask_file)

    # Add versioning information
    cice.mask_f.inputfile = f"{ocean_mask}"
    cice.mask_f.inputfile_md5 = md5sum(ocean_mask)
    cice.mask_f.history_command = f"python make_CICE_grid.py {ocean_hgrid} {ocean_mask}"

    # Add the typical crs (i.e. WGS84/EPSG4326 , but in radians).
    crs = cice.mask_f.createVariable("crs", "S1")
    crs.grid_mapping_name = "tripolar_latitude_longitude"
    crs.crs_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["radians",1,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

    cice.write_mask()


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("ocean_hgrid", help="ocean_hgrid.nc file")
    parser.add_argument("ocean_mask", help="ocean_mask.nc file")
    # to-do: add argument for CRS & output filenames?

    args = vars(parser.parse_args())

    sys.exit(cice_from_mom(**args))
