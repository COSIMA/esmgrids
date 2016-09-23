
class CICEGrid:
    """
    Make a new CICE grid and land sea mask using the MOM5 ocean grid definition
    and land sea mask.
    """

    def __init__(self, mom_grid, output_dir):
        """
        All variables are calculated here. 
        """

        self.output_dir = output_dir
        self.mom_grid = mom_grid
        self.mask_filename = os.path.join(self.output_dir, 'kmt.nc')
        self.grid_filename = os.path.join(self.output_dir, 'grid.nc')

        # Copy some stuff from MOM. FIXME: better/different way to do this.
        # Perhaps calculate here instead of in mom. Inherit from MOM? 
        self.area_t = mom_grid.area_t
        self.area_u = mom_grid.area_u

        self.x_u = mom_grid.x_u
        self.y_u = mom_grid.y_u
        self.x_t = mom_grid.x_t
        self.y_t = mom_grid.y_t

        self.clat = mom_grid.clat
        self.clon = mom_grid.clon

        self.num_lon_points = self.x_t.shape[1]
        self.num_lat_points = self.x_t.shape[0]

        # The ocean grid is double density. Containing T and U points.  The
        # cice grid is single density with separate variables for the T and
        # U points. 
        self.htn = mom_grid.dx[2::2, 0::2] + mom_grid.dx[2::2, 1::2]
        self.hte = mom_grid.dy[0::2, 2::2] + mom_grid.dy[1::2, 2::2]

        self.angle = mom_grid.angle_dx[2::2,0:-1:2]
        self.angleT = mom_grid.angle_dx[1::2,1::2]

        with nc.Dataset(self.mom_grid.mask_filename) as f:
            self.mask = np.copy(f.variables['mask'])


    def write_grid(self):
        """
        __init__() has done all the work of calculating fields. Here we make
        the netcdf file and copy over. 
        """

        f_ice = nc.Dataset(self.grid_filename, 'w')

        # Create some dimensions. 
        f_ice.createDimension('nx', self.num_lon_points)
        f_ice.createDimension('ny', self.num_lat_points)
        f_ice.createDimension('nc', 4)

        # Make all CICE grid variables. 
        ulat = f_ice.createVariable('ulat', 'f8', dimensions=('ny', 'nx'))
        ulat.units = "radians"
        ulat.title = "Latitude of U points"
        ulon = f_ice.createVariable('ulon', 'f8', dimensions=('ny', 'nx'))
        ulon.units = "radians"
        ulon.title = "Longitude of U points"
        tlat = f_ice.createVariable('tlat', 'f8', dimensions=('ny', 'nx'))
        tlat.units = "radians"
        tlat.title = "Latitude of T points"
        tlon = f_ice.createVariable('tlon', 'f8', dimensions=('ny', 'nx'))
        tlon.units = "radians"
        tlon.title = "Longitude of T points"
        htn = f_ice.createVariable('htn', 'f8', dimensions=('ny', 'nx'))
        htn.units = "cm"
        htn.title = "Width of T cells on North side."
        hte = f_ice.createVariable('hte', 'f8', dimensions=('ny', 'nx'))
        hte.units = "cm"
        hte.title = "Width of T cells on East side."
        angle = f_ice.createVariable('angle', 'f8', dimensions=('ny', 'nx'))
        angle.units = "radians"
        angle.title = "Rotation angle of U cells."
        angleT = f_ice.createVariable('angleT', 'f8', dimensions=('ny', 'nx'))
        angleT.units = "radians"
        angleT.title = "Rotation angle of T cells."
        area_t = f_ice.createVariable('tarea', 'f8', dimensions=('ny', 'nx'))
        area_t.units = "m^2"
        area_t.title = "Area of T cells."
        area_u = f_ice.createVariable('uarea', 'f8', dimensions=('ny', 'nx'))
        area_u.units = "m^2"
        area_u.title = "Area of U cells."

        area_t[:] = self.area_t[:]
        area_u[:] = self.area_u[:]

        # Now convert units: degrees -> radians. 
        tlat[:] = np.deg2rad(self.y_t)
        tlon[:] = np.deg2rad(self.x_t)
        ulat[:] = np.deg2rad(self.y_u)
        ulon[:] = np.deg2rad(self.x_u)

        # Convert from m to cm. 
        htn[:] = self.htn * 100.
        hte[:] = self.hte * 100.

        angle[:] = np.deg2rad(self.angle[:])
        angleT[:] = np.deg2rad(self.angleT[:])

        f_ice.close()


    def write_mask(self):

        input = self.mom_grid.mask_filename

        # ncrename will update the file history. 
        shutil.copyfile(input, self.mask_filename)
        with nc.Dataset(self.mask_filename, 'r+') as f:
            f.renameVariable('mask', 'kmt')

    def write(self):
        self.write_mask()
        self.write_grid()


