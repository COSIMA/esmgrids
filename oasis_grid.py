
def normalise_lons(lons):
    """
    Make a copy of lons that is between -180 and 180.
    """

    lons_n = np.copy(lons)

    lons_n[lons_n > 180] = lons_n[lons_n > 180] - 360
    lons_n[lons_n < -180] = lons_n[lons_n < -180] + 360
    return lons_n


def oasis_to_2d_corners(input_clons, input_clats):
    """
    Change from an oasis corners convention (3D) to 2D convention. Also
    changes longitudes to be from -180 to 180.

    The input_array has shape (4, rows, columns) where the first dimension is
    the corner number with 0 being lower left and increasing in
    counter-clockwise dimension:

    3---2
    |   |
    0---1

    From this see that the top corners of one grid cell (3, 2) will be the
    bottom corners (0, 1) of the grid cell above. The conversion will not
    work if this assumption does not hold.

    3---2 3---2
    |   | |   |
    0---1 0---1
    3---2 3---2
    |   | |   |
    0---1 0---1

    The output array has shape (rows + 1, columns). It contains all the 0
    corners as well as the top row 3 corners. An extra column is not needed
    because it is periodic, so the the right hand corners and the right-most
    cell are in fact the left-most corners. 
    """

    def check_array(array):
        """
        Ensure that input array supports assumption that there are no 'gaps'
        between cells. These are the assumptions explained above. 
        """

        # 0 and 3 corners have the same longitude along rows.
        assert(np.sum(array[0, 1:, :] == array[3, 0:-1, :]) ==
               array[0, 1:, :].size)
        # 1 and 2 corners have the same longitude along rows.
        assert(np.sum(array[1, 1:, :] == array[2, 0:-1, :]) ==
               array[1, 1:, :].size)
        # 0 and 1 corners have the same latitude along columns.
        assert(np.sum(array[0, :, 1:] == array[1, :, 0:-1]) ==
               array[0, :, 1:].size)
        # 2 and 3 corners have the same latitude along columns.
        assert(np.sum(array[3, :, 1:] == array[2, :, 0:-1]) ==
               array[3, :, 1:].size)
        # Wraps around in the longitude direction.
        n_array = (array + 360) % 360
        assert(np.sum(n_array[0, :, 0] == n_array[1, :, -1] ) ==
               n_array[0, :, 0].size)
        assert(np.sum(n_array[3, :, 0] == n_array[2, :, -1] ) ==
               n_array[3, :, 0].size)


    input_clons = normalise_lons(input_clons)
    check_array(input_clons)
    check_array(input_clats)

    shape = (input_clons.shape[1] + 1, input_clons.shape[2])

    output_clons = np.zeros(shape)
    output_clons[:] = np.NAN
    output_clons[:-1,:] = input_clons[0,:,:] 
    output_clons[-1,:] = input_clons[3,-1,:]
    assert(np.sum(output_clons) != np.NAN)
    # All the numbers in input should also be in output. 
    # Weird that the rounding is necessary...
    assert(set(list(np.around(input_clons.flatten(), decimals=10))) ==
           set(list(np.around(output_clons.flatten(), decimals=10))))

    output_clats = np.zeros(shape)
    output_clats[:] = np.NAN
    output_clats[:-1,:] = input_clats[0,:,:] 
    output_clats[-1,:] = input_clats[3,-1,:]
    assert(np.sum(output_clats) != np.NAN)
    assert(set(list(np.around(input_clats.flatten(), decimals=10))) ==
           set(list(np.around(output_clats.flatten(), decimals=10))))


    return output_clons, output_clats 


class OASISGrid:
    """
    This class creates 3 files: 
        - areas.nc grid cell areas for atmos (t, u, v) and ice. 
        - masks.nc land sea mask for atmos (t, u, v) and ice. 
        - grids.nc lat, lon, rotation angle and corners for atmos (t, u, v) and
          ice.
    """

    def __init__(self, atm, ice, output_dir):
        self.atm = atm
        self.ice = ice

        self.areas_filename = os.path.join(output_dir, 'areas.nc')
        self.masks_filename = os.path.join(output_dir, 'masks.nc')
        self.grids_filename = os.path.join(output_dir, 'grids.nc')

    def make_areas(self):
        """
        Make netcdf file areas.nc with cice.srf, um1t.srf, um1u.srf, um1v.srf
        """

        ice = self.ice
        atm = self.atm

        f = nc.Dataset(self.areas_filename, 'w')

        f.createDimension('nyi', ice.num_lat_points)
        f.createDimension('nxi', ice.num_lon_points)
        f.createDimension('nyat', atm.num_lat_points)
        f.createDimension('nxat', atm.num_lon_points)
        f.createDimension('nyau', atm.num_lat_points)
        f.createDimension('nxau', atm.num_lon_points)
        f.createDimension('nyav', atm.num_lat_points - 1)
        f.createDimension('nxav', atm.num_lon_points)

        cice_srf = f.createVariable('cice.srf', 'f8', dimensions=('nyi', 'nxi'))
        cice_srf.units = 'm^2'
        cice_srf.title = 'cice grid T-cell area.'
        um1t_srf = f.createVariable('um1t.srf', 'f8',
                                    dimensions=('nyat', 'nxat'))
        um1t_srf.units = 'm^2'
        um1t_srf.title = 'um1t grid area.'
        um1u_srf = f.createVariable('um1u.srf', 'f8',
                                    dimensions=('nyau', 'nxau'))
        um1u_srf.units = 'm^2'
        um1u_srf.title = 'um1u grid area.'
        um1v_srf = f.createVariable('um1v.srf', 'f8',
                                    dimensions=('nyav', 'nxav'))
        um1v_srf.units = 'm^2'
        um1v_srf.title = 'um1v grid area.'

        cice_srf[:] = self.ice.area_t[:]
        um1t_srf[:] = atm.area_t[:]
        um1u_srf[:] = atm.area_u[:]
        um1v_srf[:] = atm.area_v[:]

        f.close()


    def make_masks(self):

        f = nc.Dataset(self.masks_filename, 'w')

        f.createDimension('ny0', self.ice.num_lat_points)
        f.createDimension('nx0', self.ice.num_lon_points)
        f.createDimension('ny1', self.atm.num_lat_points)
        f.createDimension('nx1', self.atm.num_lon_points)
        f.createDimension('ny2', self.atm.num_lat_points - 1)
        f.createDimension('nx2', self.atm.num_lon_points)

        # Make the ice mask.
        mask = f.createVariable('cice.msk', 'int32', dimensions=('ny0', 'nx0'))
        mask.units = '0/1:o/l'
        mask.title = 'Ice grid T-cell land-sea mask.'
        # Flip the mask. OASIS uses 1 = masked, 0 = unmasked.
        mask[:] = (1 - self.ice.mask[:]) 

        # Atm t mask.
        mask_t = f.createVariable('um1t.msk', 'int32', dimensions=('ny1', 'nx1'))
        mask_t.units = '0/1:o/l'
        mask_t.title = 'Atm grid T-cell land-sea mask.'
        # Build the mask using the atm land fraction. 
        mask_t[:] = np.copy(self.atm.landfrac)
        mask_t[np.where(self.atm.landfrac[:] != 1)] = 0

        # Atm u mask.
        mask_u = f.createVariable('um1u.msk', 'int32', dimensions=('ny1', 'nx1'))
        mask_u.units = '0/1:o/l'
        mask_u.title = 'Atm grid U-cell land-sea mask.'
        mask_u[:] = mask_t[:]
        # Hack? make u mask by adding land points onto Western land bounary.
        for i in range(mask_u.shape[0]):
            for j in range(mask_u.shape[1] - 1):
                if mask_u[i, j+1] == 1:
                    mask_u[i, j] = 1

        # Atm v mask.
        mask_v = f.createVariable('um1v.msk', 'int32', dimensions=('ny2', 'nx2'))
        mask_v.units = '0/1:o/l'
        mask_v.title = 'Atm grid V-cell land-sea mask.'
        mask_v[:] = mask_t[:-1,:]
        # Hack? make v mask by adding land points onto Southern land bounary.
        for i in range(mask_v.shape[0] - 1):
            for j in range(mask_v.shape[1]):
                if mask_v[i+1, j] == 1:
                    mask_v[i, j] = 1

        f.close()

    def make_grids(self):
        """
        lat, lon, and corners for atmos (t, u, v) and ice.
        """

        f = nc.Dataset(self.grids_filename, 'w')

        # Ice 
        f.createDimension('nyi', self.ice.num_lat_points)
        f.createDimension('nxi', self.ice.num_lon_points)
        f.createDimension('nci', 4)
        cice_lat = f.createVariable('cice.lat', 'f8', ('nyi', 'nxi'))
        cice_lat.units = "degrees_north"
        cice_lat.title = "cice grid T-cell latitude."
        cice_lat[:] = self.ice.y_t[:]

        cice_lon = f.createVariable('cice.lon', 'f8', ('nyi', 'nxi'))
        cice_lon.units = "degrees_east"
        cice_lon.title = "cice grid T-cell longitude."
        cice_lon[:] = self.ice.x_t[:]

        cice_cla = f.createVariable('cice.cla', 'f8', ('nci', 'nyi', 'nxi'))
        cice_cla.units = "degrees_north"
        cice_cla.title = "cice grid T-cell corner latitude"
        cice_cla[:] = self.ice.clat[:]

        cice_clo = f.createVariable('cice.clo', 'f8', ('nci', 'nyi', 'nxi'))
        cice_clo.units = "degrees_east"
        cice_clo.title = "cice grid T-cell corner longitude"
        cice_clo[:] = self.ice.clon[:]

        # Atm 
        f.createDimension('nya1', self.atm.num_lat_points)
        f.createDimension('nxa1', self.atm.num_lon_points)
        f.createDimension('nya2', self.atm.num_lat_points - 1)
        f.createDimension('nxa2', self.atm.num_lon_points)
        f.createDimension('nca', 4)

        # T cells. 
        um1t_lat = f.createVariable('um1t.lat', 'f8', ('nya1', 'nxa1'))
        um1t_lat.units = "degrees_north"
        um1t_lat.title = "um1t grid center latitude"
        um1t_lat[:] = self.atm.y_t[:]

        um1t_lon = f.createVariable('um1t.lon', 'f8', ('nya1', 'nxa1'))
        um1t_lon.units = "degrees_east"
        um1t_lon.title = "um1t grid center longitude"
        um1t_lon[:] = self.atm.x_t[:]

        um1t_clat = f.createVariable('um1t.cla', 'f8', ('nca', 'nya1', 'nxa1'))
        um1t_clat.units = "degrees_north"
        um1t_clat.title = "um1t grid corner latitude"
        um1t_clat[:] = self.atm.clat_t[:]

        um1t_clon = f.createVariable('um1t.clo', 'f8', ('nca', 'nya1', 'nxa1'))
        um1t_clon.units = "degrees_east"
        um1t_clon.title = "um1t grid corner longitude"
        um1t_clon[:] = self.atm.clon_t[:]

        # U cells
        um1u_lat = f.createVariable('um1u.lat', 'f8', ('nya1', 'nxa1'))
        um1u_lat.units = "degrees_north"
        um1u_lat.title = "um1u grid center latitude"
        um1u_lat[:] = self.atm.y_u[:]

        um1u_lon = f.createVariable('um1u.lon', 'f8', ('nya1', 'nxa1'))
        um1u_lon.units = "degrees_east"
        um1u_lon.title = "um1u grid center longitude"
        um1u_lon[:] = self.atm.x_u[:]

        um1u_clat = f.createVariable('um1u.cla', 'f8', ('nca', 'nya1', 'nxa1'))
        um1u_clat.units = "degrees_north"
        um1u_clat.title = "um1u grid corner latitude"
        um1u_clat[:] = self.atm.clat_u[:]

        um1u_clon = f.createVariable('um1u.clo', 'f8', ('nca', 'nya1', 'nxa1'))
        um1u_clon.units = "degrees_east"
        um1u_clon.title = "um1u grid corner longitude"
        um1u_clon[:] = self.atm.clon_u[:]

        # V cells.
        um1v_lat = f.createVariable('um1v.lat', 'f8', ('nya2', 'nxa2'))
        um1v_lat.units = "degrees_north"
        um1v_lat.title = "um1v grid center latitude"
        um1v_lat[:] = self.atm.y_v[:]

        um1v_lon = f.createVariable('um1v.lon', 'f8', ('nya2', 'nxa2'))
        um1v_lon.units = "degrees_east"
        um1v_lon.title = "um1v grid center longitude"
        um1v_lon[:] = self.atm.x_v[:]

        um1v_clat = f.createVariable('um1v.cla', 'f8', ('nca', 'nya2', 'nxa2'))
        um1v_clat.units = "degrees_north"
        um1v_clat.title = "um1v grid corner latitude"
        um1v_clat[:] = self.atm.clat_v[:]

        um1v_clon = f.createVariable('um1v.clo', 'f8', ('nca', 'nya2', 'nxa2'))
        um1v_clon.units = "degrees_east"
        um1v_clon.title = "um1v grid corner longitude"
        um1v_clon[:] = self.atm.clon_v[:]

        f.close()


    def write(self):
        
        self.make_areas()
        self.make_masks()
        self.make_grids()


