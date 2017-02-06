
import pytest
import sh
import numpy as np
import netCDF4 as nc

data_tarball = 'test_data.tar.gz'
data_tarball_url = 'http://s3-ap-southeast-2.amazonaws.com/dp-drop/oasis-grids/test/test_data.tar.gz'

class Test():
    test_dir = os.path.dirname(os.path.realpath(__file__))
    test_data_dir = os.path.join(test_dir, 'test_data')
    test_data_tarball = os.path.join(test_dir, data_tarball)
    output_dir = os.path.join(test_data_dir, 'output')

    @pytest.fixture
    def input_dir(self):

        if not os.path.exists(self.test_data_dir):
            if not os.path.exists(self.test_data_tarball):
                sh.wget('-P', self.test_dir, data_tarball_url)
            sh.tar('zxvf', self.test_data_tarball, '-C', self.test_dir)

        return os.path.join(self.test_data_dir, 'input')

    def test_convert_mom_to_cice(self, input_dir):
        """
        Read in a MOM grid and write out a cice grid at the same resolution.
        """

        mom_grid = mom_grid.MOMGrid(os.path.join(input_dir, 'ocean_grid.nc'),

