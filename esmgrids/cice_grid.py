import numpy as np
import netCDF4 as nc

from esmgrids.base_grid import BaseGrid


class CiceGrid(BaseGrid):

    def __init__(self, **kwargs):
        self.type = "Arakawa B / C"
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

    def _create_2d_nc_var(self, f, name):
        return f.createVariable(
            name,
            "f8",
            dimensions=("ny", "nx"),
            compression="zlib",
            complevel=1,
        )

    def write(self, grid_filename, mask_filename, metadata=None):
        """
        Write out CICE grid to netcdf

        Parameters
        ----------
        grid_filename: str
            The name of the grid file to write
        mask_filename: str
            The name of the mask file to write
        metadata: dict
            Any global or variable metadata attributes to add to the files being written
        """

        # Grid file
        f = nc.Dataset(grid_filename, "w")

        # Create dimensions.
        f.createDimension("nx", self.num_lon_points)
        # nx is the grid_longitude but doesn't have a value other than its index
        f.createDimension("ny", self.num_lat_points)
        # ny is the grid_latitude but doesn't have a value other than its index

        # Make all CICE grid variables.
        # names are based on https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
        f.Conventions = "CF-1.6"

        ulat = self._create_2d_nc_var(f, "ulat")
        ulat.units = "radians"
        ulat.long_name = "Latitude of U points"
        ulat.standard_name = "latitude"
        ulon = self._create_2d_nc_var(f, "ulon")
        ulon.units = "radians"
        ulon.long_name = "Longitude of U points"
        ulon.standard_name = "longitude"
        tlat = self._create_2d_nc_var(f, "tlat")
        tlat.units = "radians"
        tlat.long_name = "Latitude of T points"
        tlat.standard_name = "latitude"
        tlon = self._create_2d_nc_var(f, "tlon")
        tlon.units = "radians"
        tlon.long_name = "Longitude of T points"
        tlon.standard_name = "longitude"

        htn = self._create_2d_nc_var(f, "htn")
        htn.units = "cm"
        htn.long_name = "Width of T cells on North side."
        htn.coordinates = "ulat tlon"
        htn.grid_mapping = "crs"
        hte = self._create_2d_nc_var(f, "hte")
        hte.units = "cm"
        hte.long_name = "Width of T cells on East side."
        hte.coordinates = "tlat ulon"
        hte.grid_mapping = "crs"

        angle = self._create_2d_nc_var(f, "angle")
        angle.units = "radians"
        angle.long_name = "Rotation angle of U cells."
        angle.standard_name = "angle_of_rotation_from_east_to_x"
        angle.coordinates = "ulat ulon"
        angle.grid_mapping = "crs"
        angleT = self._create_2d_nc_var(f, "angleT")
        angleT.units = "radians"
        angleT.long_name = "Rotation angle of T cells."
        angleT.standard_name = "angle_of_rotation_from_east_to_x"
        angleT.coordinates = "tlat tlon"
        angleT.grid_mapping = "crs"

        area_t = self._create_2d_nc_var(f, "tarea")
        area_t.units = "m^2"
        area_t.long_name = "Area of T cells."
        area_t.standard_name = "cell_area"
        area_t.coordinates = "tlat tlon"
        area_t.grid_mapping = "crs"
        area_u = self._create_2d_nc_var(f, "uarea")
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

        # Add the typical crs (i.e. WGS84/EPSG4326 , but in radians).
        crs = f.createVariable("crs", "S1")
        crs.crs_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["radians",1,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

        # Add global metadata
        if metadata:
            for attr, meta in metadata.items():
                f.setncattr(attr, meta)

        f.close()

        # Mask file
        f = nc.Dataset(mask_filename, "w")

        f.createDimension("nx", self.num_lon_points)
        f.createDimension("ny", self.num_lat_points)
        mask = self._create_2d_nc_var(f, "kmt")
        mask.grid_mapping = "crs"
        mask.standard_name = "sea_binary_mask"

        # CICE uses 0 as masked, whereas internally we use 1 as masked.
        mask[:] = 1 - self.mask_t

        # Add the typical crs (i.e. WGS84/EPSG4326 , but in radians).
        crs = f.createVariable("crs", "S1")
        crs.crs_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["radians",1,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

        # Add global metadata
        if metadata:
            for attr, meta in metadata.items():
                f.setncattr(attr, meta)

        f.close()
