import pytest
import xarray as xr
from numpy.testing import assert_allclose
from numpy import deg2rad
from subprocess import run

from esmgrids.cice_grid import cice_from_mom

# ----------------
# test data:


class MomGridFixture:
    """Generate a sample tripole grid to use as test data"""

    def __init__(self, tmp_path):
        self.path = str(tmp_path) + "/ocean_hgrid.nc"
        self.mask_path = str(tmp_path) + "/ocean_mask.nc"

        # generate a tripolar grid as test data
        run(
            [
                "ocean_grid_generator.py",
                "-r",
                "0.25",  # 4 degree grid
                "--no_south_cap",
                "--ensure_nj_even",
                "-f",
                self.path,
            ]
        )

        # open ocean_hgrid.nc
        self.ds = xr.open_dataset(self.path)

        # an ocean mask with a arbitrary mask
        self.mask_ds = xr.Dataset()
        self.mask_ds["mask"] = (self.ds.area.coarsen(ny=2).sum().coarsen(nx=2).sum()) > 5e9
        self.mask_ds.to_netcdf(self.mask_path)


class CiceGridFixture:
    """Make the CICE grid, using script under test"""

    def __init__(self, mom_grid, tmp_path):
        self.path = str(tmp_path) + "/grid.nc"
        self.kmt_path = str(tmp_path) + "/kmt.nc"
        cice_from_mom(mom_grid.path, mom_grid.mask_path,self.path, self.kmt_path)
        self.ds = xr.open_dataset(self.path, decode_cf=False)
        self.kmt_ds = xr.open_dataset(self.kmt_path, decode_cf=False)


# pytest doesn't support class fixtures, so we need these two constructor funcs
@pytest.fixture
def mom_grid(tmp_path):
    return MomGridFixture(tmp_path)


@pytest.fixture
def cice_grid(mom_grid, tmp_path):
    return CiceGridFixture(mom_grid, tmp_path)


@pytest.fixture
def test_grid_ds(mom_grid):
    # this generates the expected answers
    # In simple terms the MOM supergrid has four cells for each model grid cell. The MOM supergrid includes all edges (left and right) but CICE only uses right/east edges. (e.g. For points/edges of first cell: 0,0 is SW corner, 1,1 is middle of cell, 2,2, is NE corner/edges)

    ds = mom_grid.ds

    # u points are at top-right (NE) corner
    u_points = ds.isel(nxp=slice(2, None, 2), nyp=slice(2, None, 2))

    # t points are in middle of cell
    t_points = ds.isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))

    test_grid = xr.Dataset()

    test_grid["ulat"] = deg2rad(u_points.y)
    test_grid["ulon"] = deg2rad(u_points.x)
    test_grid["tlat"] = deg2rad(t_points.y)
    test_grid["tlon"] = deg2rad(t_points.x)

    test_grid["angle"] = deg2rad(u_points.angle_dx)  # angle at u point
    test_grid["angleT"] = deg2rad(t_points.angle_dx)

    # length of top (northern) edge of cells
    test_grid["htn"] = ds.dx.isel(nyp=slice(2, None, 2)).coarsen(nx=2).sum() * 100
    # length of right (eastern) edge of cells
    test_grid["hte"] = ds.dy.isel(nxp=slice(2, None, 2)).coarsen(ny=2).sum() * 100

    # area of cells
    test_grid["tarea"] = ds.area.coarsen(ny=2).sum().coarsen(nx=2).sum()

    # uarea is area of a cell centred around the u point
    # we need to fold on longitude and wrap in latitude to calculate this
    # drop the bottom row, new top row is reverse of current top row
    area_folded = xr.concat([ds.area.isel(ny=slice(1, None)), ds.area.isel(ny=-1, nx=slice(-1, None, -1))], dim="ny")

    # drop the first column, make the new last column the first column
    area_wrapped = xr.concat([area_folded.isel(nx=slice(1, None)), area_folded.isel(nx=0)], dim="nx")

    test_grid["uarea"] = area_wrapped.coarsen(ny=2).sum().coarsen(nx=2).sum()

    return test_grid


# ----------------
# the tests in earnest:


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_cice_var_list(cice_grid, test_grid_ds):
    # Test : Are there missing vars in cice_grid?
    assert set(test_grid_ds.variables).difference(cice_grid.ds.variables) == set()


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_cice_grid(cice_grid, test_grid_ds):
    # Test : Is the data the same as the test_grid
    for jVar in test_grid_ds.variables:
        assert_allclose(cice_grid.ds[jVar], test_grid_ds[jVar], rtol=1e-13, verbose=True, err_msg=f"{jVar} mismatch")


def test_cice_kmt(mom_grid, cice_grid):
    # Test : does the mask match
    mask = mom_grid.mask_ds.mask
    kmt = cice_grid.kmt_ds.kmt

    assert_allclose(mask, kmt, rtol=1e-13, verbose=True, err_msg="mask mismatch")


def test_cice_grid_attributes(cice_grid):
    # Test: do the expected attributes to exist in the cice ds
    cf_attributes = {
        "ulat": {"standard_name": "latitude", "units": "radians"},
        "ulon": {"standard_name": "longitude", "units": "radians"},
        "tlat": {"standard_name": "latitude", "units": "radians"},
        "tlon": {"standard_name": "longitude", "units": "radians"},
        "uarea": {
            "standard_name": "cell_area",
            "units": "m^2",
            "grid_mapping": "crs",
            "coordinates": "ulat ulon",
        },
        "tarea": {
            "standard_name": "cell_area",
            "units": "m^2",
            "grid_mapping": "crs",
            "coordinates": "tlat tlon",
        },
        "angle": {
            "standard_name": "angle_of_rotation_from_east_to_x",
            "units": "radians",
            "grid_mapping": "crs",
            "coordinates": "ulat ulon",
        },
        "angleT": {
            "standard_name": "angle_of_rotation_from_east_to_x",
            "units": "radians",
            "grid_mapping": "crs",
            "coordinates": "tlat tlon",
        },
        "htn": {"units": "cm", "coordinates": "ulat tlon", "grid_mapping": "crs"},
        "hte": {"units": "cm", "coordinates": "tlat ulon", "grid_mapping": "crs"},
    }

    for iVar in cf_attributes.keys():
        print(cice_grid.ds[iVar])

        for jAttr in cf_attributes[iVar].keys():
            assert cice_grid.ds[iVar].attrs[jAttr] == cf_attributes[iVar][jAttr]


def test_crs_exist(cice_grid):
    # Test: has the crs been added ?
    # todo: open with GDAL and rioxarray and confirm they find the crs?
    assert hasattr(cice_grid.ds, "crs")
    assert hasattr(cice_grid.kmt_ds, "crs")


def test_inputs_logged(cice_grid, mom_grid):
    # Test: have the source data been logged ?

    for ds in [cice_grid.ds, cice_grid.kmt_ds]:
        assert hasattr(ds, "inputfile"), "inputfile attribute missing"
        assert hasattr(ds, "inputfile_md5"), "inputfile md5sum attribute missing"

        sys_md5 = run(["md5sum", ds.inputfile], capture_output=True, text=True)
        sys_md5 = sys_md5.stdout.split(" ")[0]
        assert ds.inputfile_md5 == sys_md5, f"inputfile md5sum attribute incorrect, {ds.inputfile_md5} != {sys_md5}"

    assert (
        cice_grid.ds.inputfile == mom_grid.path
    ), "inputfile attribute incorrect ({cice_grid.ds.inputfile} != {mom_grid.path})"
    assert (
        cice_grid.kmt_ds.inputfile == mom_grid.mask_path
    ), "mask inputfile attribute incorrect ({cice_grid.kmt_ds.inputfile} != {mom_grid.mask_path})"
