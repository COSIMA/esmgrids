
import pytest
import sh
import os
import numpy as np
import netCDF4 as nc

from esmgrids.mom_grid import MomGrid
from esmgrids.cice_grid import CiceGrid

data_tarball = 'test_data.tar.gz'
data_tarball_url = 'http://s3-ap-southeast-2.amazonaws.com/dp-drop/esmgrids/test/test_data.tar.gz'

class Test():
    test_dir = os.path.dirname(os.path.realpath(__file__))
    test_data_dir = os.path.join(test_dir, 'test_data')
    test_data_tarball = os.path.join(test_dir, data_tarball)
    out_dir = os.path.join(test_data_dir, 'output')

    @pytest.fixture
    def output_dir(self):
        return self.out_dir

    @pytest.fixture
    def input_dir(self):

        if not os.path.exists(self.test_data_dir):
            if not os.path.exists(self.test_data_tarball):
                sh.wget('-P', self.test_dir, data_tarball_url)
            sh.tar('zxvf', self.test_data_tarball, '-C', self.test_dir)

        return os.path.join(self.test_data_dir, 'input')

    def test_convert_mom_to_cice(self, input_dir, output_dir):
        """
        Read in a MOM grid and write out a cice grid at the same resolution.
        """

        mask = os.path.join(input_dir, 'ocean_01_mask.nc')
        hgrid = os.path.join(input_dir, 'ocean_01_hgrid.nc')
        mom = MomGrid.fromfile(hgrid, mask_file=mask)
        cice = CiceGrid.fromgrid(mom)
        grid_file = os.path.join(output_dir, 'cice_grid.nc')
        mask_file = os.path.join(output_dir, 'cice_mask.nc')
        cice.write(grid_file, mask_file)

        # FIXME tests for the CICE grid.
