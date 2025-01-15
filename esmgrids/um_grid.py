class UMGrid:
    """
    Interpolate the ocean mask onto the UM grid.

    Creates the land fraction and mask, later used as an input to the UM.
    """

    def __init__(self, um_restart, num_lon_points, num_lat_points, mom_grid, output_dir):

        self.SOUTHERN_EXTENT = -89.995002746582031
        self.SOUTHERN_EXTENT_CNR = -89.999496459960938
        self.NORTHERN_EXTENT = 89.995002746582031
        self.NORTHERN_EXTENT_CNR = 89.999496459960938

        self.mom_grid = mom_grid
        self.um_restart = um_restart

        self.um_restart_output = os.path.join(output_dir, os.path.basename(um_restart))
        self.lfrac_filename_nc = os.path.join(output_dir, "lfrac.nc")
        self.lfrac_filename_um = os.path.join(output_dir, "lfrac")
        self.mask_filename_nc = os.path.join(output_dir, "qrparm.mask.nc")
        self.mask_filename_um = os.path.join(output_dir, "qrparm.mask")

        self.num_lon_points = num_lon_points
        self.num_lat_points = num_lat_points
        self.corners = 4

        # Set lats and lons.
        self.lon = np.linspace(0, 360, num_lon_points, endpoint=False)
        self.lat = np.linspace(-90, 90, num_lat_points)
        dx_half = 360.0 / num_lon_points / 2.0
        dy_half = 180.0 / (num_lat_points - 1) / 2.0

        # Similar to lon, lat but specify the coordinate at every grid
        # point. Also it wraps along longitude.
        self.x_t = np.tile(self.lon, (num_lat_points, 1))
        self.y_t = np.tile(self.lat, (num_lon_points, 1))
        self.y_t = self.y_t.transpose()

        self.x_u = self.x_t + dx_half
        self.y_u = self.y_t
        self.x_v = self.x_t
        self.y_v = self.y_t + dy_half

        def make_corners(x, y, dx, dy):

            # Set grid corners, we do these one corner at a time. Start at the
            # bottom left and go anti-clockwise. This is the OASIS convention.
            clon = np.empty((self.corners, x.shape[0], x.shape[1]))
            clon[:] = np.nan
            clon[0, :, :] = x - dx
            clon[1, :, :] = x + dx
            clon[2, :, :] = x + dx
            clon[3, :, :] = x - dx
            assert not np.isnan(np.sum(clon))

            clat = np.empty((self.corners, x.shape[0], x.shape[1]))
            clat[:] = np.nan
            clat[0, :, :] = y[:, :] - dy
            clat[1, :, :] = y[:, :] - dy
            clat[2, :, :] = y[:, :] + dy
            clat[3, :, :] = y[:, :] + dy

            # The bottom latitude band should always be Southern extent, for
            # all t, u, v.
            clat[0, 0, :] = -90
            clat[1, 0, :] = -90

            # The top latitude band should always be Northern extent, for all
            # t, u, v.
            clat[2, -1, :] = 90
            clat[3, -1, :] = 90

            assert not np.isnan(np.sum(clat))

            return clon, clat

        self.clon_t, self.clat_t = make_corners(self.x_t, self.y_t, dx_half, dy_half)
        self.clon_u, self.clat_u = make_corners(self.x_u, self.y_u, dx_half, dy_half)
        self.clon_v, self.clat_v = make_corners(self.x_v, self.y_v, dx_half, dy_half)

        # The Northerly v points are going to be beyond the domain. Remove them
        self.y_v = self.y_v[:-1, :]
        self.x_v = self.x_v[:-1, :]
        self.clat_v = self.clat_v[:, :-1, :]
        self.clon_v = self.clon_v[:, :-1, :]
        # self.area_v = self.area_v[:-1, :]

        # Now that the grid is made we fix it up. We don't go from -90 to 90
        # but from self.SOUTHERN_EXTENT to self.NORTHERN_EXTENT. As far as I
        # can tell this is due to the SCRIP remapping library not handling the
        # poles properly and making bad weights. There is a test for this in
        # tests/test_scrip_remapping.py. If the tests don't pass there's no
        # point running the model with those remapping files.
        def fix_grid():
            self.lat[0] = self.SOUTHERN_EXTENT
            self.lat[-1] = self.NORTHERN_EXTENT
            self.y_t[0, :] = self.SOUTHERN_EXTENT
            self.y_t[-1, :] = self.NORTHERN_EXTENT
            self.y_u[0, :] = self.SOUTHERN_EXTENT
            self.y_u[-1, :] = self.NORTHERN_EXTENT

            def fix_corners(clat):
                clat[0, 0, :] = self.SOUTHERN_EXTENT_CNR
                clat[1, 0, :] = self.SOUTHERN_EXTENT_CNR
                clat[2, -1, :] = self.NORTHERN_EXTENT_CNR
                clat[3, -1, :] = self.NORTHERN_EXTENT_CNR

            fix_corners(self.clat_t)
            fix_corners(self.clat_u)
            fix_corners(self.clat_v)

        fix_grid()

        # Use corners to calculate areas.
        self.area_t = self.calc_area(self.clon_t, self.clat_t)
        self.area_u = self.calc_area(self.clon_u, self.clat_u)
        self.area_v = self.calc_area(self.clon_v, self.clat_v)

        # This is defined after a call to make_landfrac.
        self.landfrac = None
        self.mask_t = None
        self.mask_u = None
        self.mask_v = None

    def calc_area(self, clons, clats):
        """
        Calculate the area of lat-lon polygons.

        We project sphere onto a flat surface using an equal area projection
        and then calculate the area of flat polygon.
        """

        def area_polygon(p):
            """
            Calculate the area of a polygon.

            Input is a polygon represented as a list of (x,y) vertex
            coordinates, implicitly wrapping around from the last vertex to the
            first.

            See http://stackoverflow.com/questions/451426/
                how-do-i-calculate-the-surface-area-of-a-2d-polygon
            """

            def segments(v):
                return zip(v, v[1:] + [v[0]])

            return 0.5 * abs(sum(x0 * y1 - x1 * y0 for ((x0, y0), (x1, y1)) in segments(p)))

        areas = np.zeros_like(clons[0])
        areas[:] = np.nan

        m = Basemap(
            projection="laea",
            resolution="h",
            llcrnrlon=0,
            llcrnrlat=-90.0,
            urcrnrlon=360,
            urcrnrlat=90.0,
            lat_0=-90,
            lon_0=0,
        )

        x, y = m(clons, clats)

        for i in range(x[0, :].shape[0]):
            for j in range(x[0, :].shape[1]):
                areas[i, j] = area_polygon(zip(x[:, i, j], y[:, i, j]))

        assert np.sum(areas) is not np.nan
        assert np.min(areas) > 0
        assert abs(1 - np.sum(areas) / EARTH_AREA) < 2e-4

        return areas

    def make_antarctic_mask(self, southern_lat, grid_lats):
        """
        Create mask on grid_lats to mask out everything South of a particular
        lat.
        """

        def find_nearest_larger(val, array):
            """
            Find the value which is nearest and larger than val in array.
            """

            s_array = np.sort(array, axis=None)
            r = np.searchsorted(s_array, val, side="right")
            return s_array[r]

        mask = np.zeros_like(grid_lats, dtype=bool)

        # Find Southern latitude of the destination that matches source.
        closest = find_nearest_larger(southern_lat, grid_lats)
        excluded_row = np.where(closest == grid_lats)[0][0]

        # Expect that lower latitudes have lower indices.
        assert all(grid_lats[excluded_row] > grid_lats[excluded_row - 1])
        # Mask out all latitude bands equal to and less than closest.
        mask[0:excluded_row, :] = True

        return mask

    def make_landfrac(self):
        """
        Regrid the ocean mask to create new land-sea fraction.
        """

        src_clons, src_clats = oasis_to_2d_corners(self.mom_grid.clon, self.mom_grid.clat)
        dest_clons, dest_clats = oasis_to_2d_corners(self.clon_t, self.clat_t)

        # The source grid is not defined South of -81. The easiest way to
        # deal with this is to mask out the destination during regridding
        # and then set it all to land.
        ant_mask = self.make_antarctic_mask(np.min(self.mom_grid.y_t), self.y_t)

        # Set up regridder with source and destination grid defs. All lons are
        # normalised -180, 180
        src_lons = normalise_lons(self.mom_grid.x_t)
        dest_lons = normalise_lons(self.x_t)

        r = Regridder(
            src_lons,
            self.mom_grid.y_t,
            src_clons,
            src_clats,
            None,
            dest_lons,
            self.y_t,
            dest_clons,
            dest_clats,
            ant_mask,
        )

        # Do regridding of mom ocean mask. This will result in an
        # 'ocean fraction' not a land fraction.
        self.landfrac = r.regrid(self.mom_grid.mask)

        # Check regridding, ensure that src and dest masses are close.
        src_mass = np.sum(self.mom_grid.area_t * self.mom_grid.mask)
        dest_mass = np.sum(self.area_t * self.landfrac)
        # assert(np.isclose(1, src_mass / dest_mass, atol=1e-5))
        # FIXME: this is not very close!
        assert np.isclose(1, src_mass / dest_mass, atol=1e-3)

        # The destination has been masked out over Antarctica for regridding
        # purposes, set that area to land.
        self.landfrac[np.where(ant_mask)] = 0

        # Flip so that we have land fraction, rather than ocean fraction.
        self.landfrac[:] = abs(1 - self.landfrac[:])
        # Clean up points which have a very small land fraction.
        self.landfrac[np.where(self.landfrac[:] < 0.01)] = 0
        self.landfrac[np.where(self.landfrac[:] > 1)] = 1

    def put_basic_header(self, file):
        """
        Put in the basic netcdf header elements: lat, lon, time.
        """

        file.createDimension("longitude", self.num_lon_points)
        file.createDimension("latitude", self.num_lat_points)
        file.createDimension("t")

        lon = file.createVariable("longitude", "f8", dimensions=("longitude"))
        lon.long_name = "longitude"
        lon.standard_name = "longitude"
        lon.units = "degrees_east"
        lon.point_spacing = "even"
        lon.module = ""

        lat = file.createVariable("latitude", "f8", dimensions=("latitude"))
        lat.long_name = "latitude"
        lat.standard_name = "latitude"
        lat.units = "degrees_north"
        lat.point_spacing = "even"

        t = file.createVariable("t", "f8", dimensions=("t"))
        t.long_name = "t"
        t.units = "days since 0001-01-01 00:00:00"
        t.time_origin = "01-JAN-0001:00:00:00"
        t[0] = 0

        lon[:] = self.lon
        lat[:] = self.lat

    def write_landfrac(self, convert_to_um=False):
        """
        Write out the land fraction.
        """

        assert self.landfrac is not None

        f = nc.Dataset(self.lfrac_filename_nc, "w", format="NETCDF3_CLASSIC")

        # Put in basic header elements lat, lon, time etc.
        self.put_basic_header(f)
        f.createDimension("ht", 1)

        ht = f.createVariable("ht", "f8", dimensions=("ht"))
        ht.long_name = "Height"
        ht.units = "m"
        ht.positive = "up"
        lsm = f.createVariable("lsm", "f8", dimensions=("t", "ht", "latitude", "longitude"))
        lsm.name = "lsm"
        lsm.title = "Stash code = 505"
        lsm.title = "Land fraction in grid box"
        lsm.valid_min = 0.0
        lsm.valid_max = 1.0

        lsm[0, 0, :, :] = self.landfrac[:]
        f.close()

        # Convert to UM format.
        if convert_to_um:
            mkancil = Mkancil()
            ret = mkancil.convert_lfrac()
            assert ret == 0
            assert os.path.exists(self.lfrac_filename_um)

    def write_mask(self, convert_to_um=False):
        """
        Write out mask used by the UM.

        This mask is used to differentiate between points that have some land
        fraction and those which have none at all.
        """
        assert self.landfrac is not None

        f = nc.Dataset(self.mask_filename_nc, "w", format="NETCDF3_CLASSIC")
        # Put in basic header elements lat, lon, time etc.
        self.put_basic_header(f)

        f.createDimension("surface", 1)

        surface = f.createVariable("surface", "f8", dimensions=("surface"))
        surface.long_name = "Surface"
        surface.units = "level"
        surface.positive = "up"

        lsm = f.createVariable("lsm", "f8", dimensions=("t", "surface", "latitude", "longitude"))
        lsm.name = "lsm"
        lsm.title = "LAND MASK (No halo) (LAND=TRUE)"
        lsm.valid_min = 0.0
        lsm.valid_max = 1.0

        # Make the mask using the land fraction.
        mask = np.copy(self.landfrac)
        mask[np.where(self.landfrac[:] != 0)] = 1
        lsm[0, 0, :, :] = mask[:]
        f.close()

        # Convert to UM format.
        if convert_to_um:
            mkancil = Mkancil()
            ret = mkancil.convert_mask()
            assert ret == 0
            assert os.path.exists(self.mask_filename_um)

    def write(self):

        self.write_landfrac()
        self.write_mask()

        shutil.copyfile(self.um_restart, self.um_restart_output)

        # Update the um restart with new mask and landfrac.
        with nc.Dataset(self.mask_filename_nc) as mask_f:
            mask = np.copy(mask_f.variables["lsm"][0, 0, :, :])
            # Flip because we use True to mean masked, UM uses True to mean
            # land.
            mask = abs(1 - mask)

            with nc.Dataset(self.lfrac_filename_nc) as lfrac_f:
                lfrac = lfrac_f.variables["lsm"][:]
                remask(self.um_restart_output, mask, lfrac)
