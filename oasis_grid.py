
from base_grid import BaseGrid

class OasisGrid(BaseGrid)
    """
    Python representation of OASIS grid including:
        - grid cell centre and corners
        - grid cell area
        - grid cell masking
    """
    def __init__(self, grid_name, model_grid, cells):

        # OASIS wants grid names to be 4 characters long and we to reserve one
        # character of the 't' or 'u'.
        assert(len(grid.name) >= 3)

        self.name = grid_name
        # Arakawa B or C grids.
        self.grid_type = model_grid.type
        assert(self.grid_type == 'Arakawa B' or self.grid_type == 'Arakawa C')

        self.model_grid = model_grid
        self.cells = cells

    def write_grids(self, grids_filename):
        """
        Make netcdf file grids.nc
        """

        f = nc.Dataset(grids_filename, 'r+')

        for cell in in self.cells:
            assert(cell == 't' or cell == 'u' or cell == 'v')

            ny_dim = 'ny{}_{}'.format(cell, self.name)
            nx_dim = 'nx{}_{}'.format(cell, self.name)
            nc_dim = 'nc{}_{}'.format(cell, self.name)
            if cell == 't':
                f.createDimension(ny_dim, self.model_grid.num_lat_points)
            else:
                f.createDimension(ny_dim, self.model_grid.num_lat_points - 1)
            f.createDimension(nx_dim, self.model_grid.num_lon_points)
            f.createDimension(nc_dim, 4)

            lat_var = '{}{}.lat'.format(self.name[:3], cell)
            lon_var = '{}{}.lon'.format(self.name[:3], cell)
            cla_var = '{}{}.cla'.format(self.name[:3], cell)
            clo_var = '{}{}.clo'.format(self.name[:3], cell)

            if cell == 't':
                y = self.model_grid.y_t[:]
                x = self.model_grid.x_t[:]
                clat = self.model_grid.clat_t[:]
                clon = self.model_grid.clon_t[:]
            elif cell == 'u':
                y = self.model_grid.y_u[:]
                x = self.model_grid.x_u[:]
                clat = self.model_grid.clat_u[:]
                clon = self.model_grid.clon_u[:]
            else:
                assert(cell == 'v')
                assert(self.grid_type == 'Arakawa C')
                y = self.model_grid.y_v[:]
                x = self.model_grid.x_v[:]
                clat = self.model_grid.clat_v[:]
                clon = self.model_grid.clon_v[:]

            tmp = f.createVariable(lat_var, 'f8', (ny_dim, nx_dim))
            tmp.units = "degrees_north"
            tmp.title = "{} grid {}-cell latitude.".format(self.name, cell.upper())
            tmp[:] = self.model_grid.y[:]

            tmp = f.createVariable(lon_var, 'f8', (ny_dim, nx_dim))
            tmp.units = "degrees_east"
            tmp.title = "{} grid {}-cell longitude.".format(self.name, cell.upper())
            tmp[:] = self.model_grid.x[:]

            tmp = f.createVariable(cla_var, 'f8', (nc_dim, ny_dim, nx_dim))
            tmp.units = "degrees_north"
            tmp.title = "{} grid {}-cell corner latitude".format(self.name, cell.upper())
            tmp[:] = self.model_grid.clat[:]

            tmp = f.createVariable(clo_var, 'f8', (nc_dim, ny_dim, nx_dim))
            tmp.units = "degrees_east"
            tmp.title = "{} grid {}-cell corner longitude".format(self.name, cell.upper())
            tmp[:] = self.model_grid.clon[:]

        f.close()


    def write_areas(self, areas_filename):
        """
        Make netcdf file areas.nc
        """

        f = nc.Dataset(self.areas_filename, 'r+')

        for cell in self.cells:
            assert(cell == 't' or cell == 'u' or cell == 'v')

            ny_dim = 'ny{}_{}'.format(cell, self.name)
            nx_dim = 'nx{}_{}'.format(cell, self.name)
            if cell == 't':
                f.createDimension(ny_dim, self.model_grid.num_lat_points)
            else:
                f.createDimension(ny_dim, self.model_grid.num_lat_points - 1)
            f.createDimension(nx_dim, self.model_grid.num_lon_points)

            srf_var = '{}{}.srf'.format(self.name[:3], cell)

            if cell == 't':
                area = self.model_grid.area_t[:]
            elif cell == 'u':
                area = self.model_grid.area_u[:]
            else:
                assert(cell == 'v')
                assert(self.grid_type == 'Arakawa C')
                area = self.model_grid.area_v[:]

            tmp = f.createVariable(lat_var, 'f8', (ny_dim, nx_dim))
            tmp.units = "m^2"
            tmp.title = "{} grid {}-cell area.".format(self.name, cell.upper())
            tmp[:] = area[:]

        f.close()


    def write_masks(self, masks_filename):

        f = nc.Dataset(self.masks_filename, 'r+')

        for cell in self.cells:
            assert(cell == 't' or cell == 'u' or cell == 'v')

            ny_dim = 'ny{}_{}'.format(cell, self.name)
            nx_dim = 'nx{}_{}'.format(cell, self.name)
            if cell == 't':
                f.createDimension(ny_dim, self.model_grid.num_lat_points)
            else:
                f.createDimension(ny_dim, self.model_grid.num_lat_points - 1)
            f.createDimension(nx_dim, self.model_grid.num_lon_points)

            msk_var = '{}{}.msk'.format(self.name[:3], cell)

            if cell == 't':
                mask = self.model_grid.mask_t[:]
            elif cell == 'u':
                mask = self.model_grid.mask_u[:]
            else:
                assert(cell == 'v')
                assert(self.grid_type == 'Arakawa C')
                mask = self.model_grid.mask_v[:]

            tmp = f.createVariable(msk_var, 'f8', (ny_dim, nx_dim))
            tmp.units = '0/1:o/l'
            tmp.title = "{} grid {}-cell land-sea mask.".format(self.name, cell.upper())
            # Flip the mask. OASIS uses 1 = masked, 0 = unmasked.
            tmp[:] = (1 - mask[:])

        f.close()
