import pytest
import sh
import os
import sys
import numpy as np

my_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(my_dir, "../"))

from esmgrids.core2_grid import Core2Grid  # noqa
from esmgrids.util import calc_area_of_polygons  # noqa

data_tarball = "test_data.tar.gz"
data_tarball_url = "http://s3-ap-southeast-2.amazonaws.com" "/dp-drop/esmgrids/test/test_data.tar.gz"

EARTH_AREA = 510072000e6


def check_for_gaps(corners, normalise_lons=False):
    """
    The top corners of one grid cell (3, 2) should be the bottom corners
    (0, 1) of the grid cell above.

    3---2 3---2
    |   | |   |
    0---1 0---1
    3---2 3---2
    |   | |   |
    0---1 0---1
    """

    # 0 and 3 corners are the same on columns
    assert (corners[0, 1:, :] == corners[3, 0:-1, :]).all()

    # 1 and 2 corners are the same on columns.
    assert (corners[1, 1:, :] == corners[2, 0:-1, :]).all()

    # 0 and 1 corners are the same on rows
    assert (corners[0, :, 1:] == corners[1, :, 0:-1]).all()

    # 2 and 3 corners are the same on rows
    assert (corners[3, :, 1:] == corners[2, :, 0:-1]).all()

    # The columns wrap around.
    n_corners = np.copy(corners)
    if normalise_lons:
        n_corners = (corners + 360) % 360

    assert (n_corners[0, :, 0] == corners[1, :, -1]).all()
    assert (n_corners[3, :, 0] == n_corners[2, :, -1]).all()


def check_corners(grid):
    """
    Do some checks on the corners.
    """

    # Check that no single cell has two corners at the same location.

    # Check total area
    area_t = calc_area_of_polygons(grid.clon_t, grid.clat_t)
    assert abs(1 - np.sum(area_t) / EARTH_AREA) < 5e-4

    # Check for gaps between cells.
    check_for_gaps(grid.clon_t, normalise_lons=True)
    check_for_gaps(grid.clat_t)


class Test:
    test_dir = os.path.dirname(os.path.realpath(__file__))
    test_data_dir = os.path.join(test_dir, "test_data")
    test_data_tarball = os.path.join(test_dir, data_tarball)
    out_dir = os.path.join(test_data_dir, "output")

    @pytest.fixture
    def output_dir(self):
        return self.out_dir

    @pytest.fixture
    def input_dir(self):

        if not os.path.exists(self.test_data_dir):
            if not os.path.exists(self.test_data_tarball):
                sh.wget("-P", self.test_dir, data_tarball_url)
            sh.tar("zxvf", self.test_data_tarball, "-C", self.test_dir)

        return os.path.join(self.test_data_dir, "input")

    @pytest.mark.broken
    def test_corners(self, input_dir):
        """
        Check some corners fields, clat_t, clon_t etc.
        """

        hgrid = os.path.join(input_dir, "t_10.0001.nc")
        core2 = Core2Grid(hgrid)
        area_t = check_corners(core2)

        # mask = os.path.join(input_dir, 'ocean_01_mask.nc')
        # hgrid = os.path.join(input_dir, 'ocean_01_hgrid.nc')
        # mom = MomGrid.fromfile(hgrid, mask_file=mask)
        # check_corners(mom)


# To be removed once we have a working test
def test_dummy():
    pass
