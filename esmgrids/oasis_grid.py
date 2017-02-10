
import os
import netCDF4 as nc
from .base_grid import BaseGrid

class OasisGrid(BaseGrid):
    """
    Python representation of OASIS grid including:
        - grid cell centre and corners
        - grid cell area
        - grid cell masking
    """
    def __init__(self, grid_name, model_grid):

        # OASIS wants grid names to be 4 characters long and we to reserve one
        # character of the 't' or 'u'.
        assert(len(grid_name) >= 3)

        self.name = grid_name
        self.grid_type = model_grid.type

        if self.grid_type == 'Arakawa B':
            self.cells = ('t', 'u')
        elif self.grid_type == 'Arakawa C':
            self.cells = ('t', 'u', 'v')
        else:
            self.cells = ('t')

        self.model_grid = model_grid

    def write_grids(self, grids_filename):
        """
        Make netcdf file grids.nc
        """

        if not os.path.exists(grids_filename):
            f = nc.Dataset(grids_filename, 'w')
        else:
            f = nc.Dataset(grids_filename, 'a')

        for cell in self.cells:
            assert(cell == 't' or cell == 'u' or cell == 'v')

            ny_dim = 'ny{}_{}'.format(cell, self.name)
            nx_dim = 'nx{}_{}'.format(cell, self.name)
            nc_dim = 'nc{}_{}'.format(cell, self.name)
            try:
                f.createDimension(ny_dim, self.model_grid.num_lat_points)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e
            try:
                f.createDimension(nx_dim, self.model_grid.num_lon_points)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e
            try:
                f.createDimension(nc_dim, 4)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e

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
                assert self.grid_type != 'Arakawa A'
                assert self.grid_type != 'Spectral'
                y = self.model_grid.y_u[:]
                x = self.model_grid.x_u[:]
                clat = self.model_grid.clat_u[:]
                clon = self.model_grid.clon_u[:]
            else:
                assert cell == 'v'
                assert self.grid_type == 'Arakawa C'
                y = self.model_grid.y_v[:]
                x = self.model_grid.x_v[:]
                clat = self.model_grid.clat_v[:]
                clon = self.model_grid.clon_v[:]

            try:
                tmp = f.createVariable(lat_var, 'f8', (ny_dim, nx_dim))
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    tmp = f.variables[lat_var]
                else:
                    raise e
            tmp.units = "degrees_north"
            tmp.title = "{} grid {}-cell latitude.".format(self.name, cell.upper())
            tmp[:] = y[:]

            try:
                tmp = f.createVariable(lon_var, 'f8', (ny_dim, nx_dim))
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    tmp = f.variables[lon_var]
                else:
                    raise e
            tmp.units = "degrees_east"
            tmp.title = "{} grid {}-cell longitude.".format(self.name, cell.upper())
            tmp[:] = x[:]

            try:
                tmp = f.createVariable(cla_var, 'f8', (nc_dim, ny_dim, nx_dim))
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    tmp = f.variables[cla_var]
                else:
                    raise e
            tmp.units = "degrees_north"
            tmp.title = "{} grid {}-cell corner latitude".format(self.name, cell.upper())
            tmp[:] = clat[:]

            try:
                tmp = f.createVariable(clo_var, 'f8', (nc_dim, ny_dim, nx_dim))
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    tmp = f.variables[clo_var]
                else:
                    raise e
            tmp.units = "degrees_east"
            tmp.title = "{} grid {}-cell corner longitude".format(self.name, cell.upper())
            tmp[:] = clon[:]

        f.close()


    def write_areas(self, areas_filename):
        """
        Make netcdf file areas.nc
        """

        if not os.path.exists(areas_filename):
            f = nc.Dataset(areas_filename, 'w')
        else:
            f = nc.Dataset(areas_filename, 'a')

        for cell in self.cells:
            assert(cell == 't' or cell == 'u' or cell == 'v')

            ny_dim = 'ny{}_{}'.format(cell, self.name)
            nx_dim = 'nx{}_{}'.format(cell, self.name)
            try:
                f.createDimension(ny_dim, self.model_grid.num_lat_points)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e
            try:
                f.createDimension(nx_dim, self.model_grid.num_lon_points)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e

            srf_var = '{}{}.srf'.format(self.name[:3], cell)

            if cell == 't':
                area = self.model_grid.area_t[:]
            elif cell == 'u':
                assert self.grid_type != 'Arakawa A'
                assert self.grid_type != 'Spectral'
                area = self.model_grid.area_u[:]
            else:
                assert(cell == 'v')
                assert(self.grid_type == 'Arakawa C')
                area = self.model_grid.area_v[:]

            try:
                tmp = f.createVariable(srf_var, 'f8', (ny_dim, nx_dim))
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    tmp = f.variables[srf_var]
                else:
                    raise e
            tmp.units = "m^2"
            tmp.title = "{} grid {}-cell area.".format(self.name, cell.upper())
            tmp[:] = area[:]

        f.close()


    def write_masks(self, masks_filename):

        if not os.path.exists(masks_filename):
            f = nc.Dataset(masks_filename, 'w')
        else:
            f = nc.Dataset(masks_filename, 'a')

        for cell in self.cells:
            assert(cell == 't' or cell == 'u' or cell == 'v')

            ny_dim = 'ny{}_{}'.format(cell, self.name)
            nx_dim = 'nx{}_{}'.format(cell, self.name)
            try:
                f.createDimension(ny_dim, self.model_grid.num_lat_points)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e
            try:
                f.createDimension(nx_dim, self.model_grid.num_lon_points)
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    pass
                else:
                    raise e

            msk_var = '{}{}.msk'.format(self.name[:3], cell)

            if cell == 't':
                mask = self.model_grid.mask_t[:]
            elif cell == 'u':
                assert self.grid_type != 'Arakawa A'
                assert self.grid_type != 'Spectral'
                mask = self.model_grid.mask_u[:]
            else:
                assert(cell == 'v')
                assert(self.grid_type == 'Arakawa C')
                mask = self.model_grid.mask_v[:]

            if len(mask.shape) == 4:
                mask = mask[0, 0, :, :]
            elif len(mask.shape) == 3:
                mask = mask[0, :, :]

            try:
                tmp = f.createVariable(msk_var, 'i4', (ny_dim, nx_dim))
            except RuntimeError as e:
                if str(e) == 'NetCDF: String match to name in use':
                    tmp = f.variables[msk_var]
                else:
                    raise e
            tmp.units = '0/1:o/l'
            tmp.title = "{} grid {}-cell land-sea mask.".format(self.name, cell.upper())
            # OASIS uses 1 = masked, 0 = unmasked. This is the same as the
            # esmgrids convention
            tmp[:] = mask[:]

        f.close()
