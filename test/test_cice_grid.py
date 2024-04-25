import pytest
import xarray as xr
from numpy.testing import assert_allclose
from numpy import deg2rad
from subprocess import run
from pathlib import Path

from esmgrids.mom_grid import MomGrid
from esmgrids.cice_grid import CiceGrid

# create test grids at 4 degrees and 0.1 degrees
# 4 degress is the lowest tested in ocean_model_grid_generator
# going higher resolution than 0.1 has too much computational cost
_test_resolutions = [4, 0.1]

# run test using the valid cice variants
_variants = ["cice5-auscom", None]


# so that our fixtures are only created once in this pytest module, we need this special version of 'tmp_path'
@pytest.fixture(scope="module")
def tmp_path(tmp_path_factory: pytest.TempdirFactory) -> Path:
    return tmp_path_factory.mktemp("temp")


# ----------------
# test data:
class MomGridFixture:
    """Generate a sample tripole grid to use as test data"""

    def __init__(self, res, tmp_path):
        self.path = str(tmp_path) + "/ocean_hgrid.nc"
        self.mask_path = str(tmp_path) + "/ocean_mask.nc"

        # generate a tripolar grid as test data
        run(
            [
                "ocean_grid_generator.py",
                "-r",
                str(1 / res),
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

    def __init__(self, mom_grid, tmp_path, variant):
        self.path = str(tmp_path) + "/grid.nc"
        self.kmt_path = str(tmp_path) + "/kmt.nc"

        run_cmd = [
            "cice_from_mom",
            "--ocean_hgrid",
            mom_grid.path,
            "--ocean_mask",
            mom_grid.mask_path,
            "--cice_grid",
            self.path,
            "--cice_kmt",
            self.kmt_path,
        ]
        if variant is not None:
            run_cmd.append("--cice_variant")
            run_cmd.append(variant)
        run(run_cmd)

        self.ds = xr.open_dataset(self.path, decode_cf=False)
        self.kmt_ds = xr.open_dataset(self.kmt_path, decode_cf=False)


def gen_grid_ds(mom_grid, variant):
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
    if variant == "cice5-auscom":
        test_grid["angleT"] = deg2rad(t_points.angle_dx)
    else:  # cice6
        test_grid["anglet"] = deg2rad(t_points.angle_dx)

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


# pytest doesn't support class fixtures, so we need these two constructor funcs
@pytest.fixture(scope="module", params=_test_resolutions)
def mom_grid(request, tmp_path):
    return MomGridFixture(request.param, tmp_path)


# the variant neews to be the same for both the cice_grid and the test_grid, so bundle them
@pytest.fixture(scope="module", params=_variants)
def grids(request, mom_grid, tmp_path):
    return {"cice": CiceGridFixture(mom_grid, tmp_path, request.param), "test_ds": gen_grid_ds(mom_grid, request.param)}


# ----------------
# the tests in earnest:


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_cice_var_list(grids):
    # Test : Are there missing vars in cice_grid?
    assert set(grids["test_ds"].variables).difference(grids["cice"].ds.variables) == set()


def test_cice_dims(grids):
    # Test : Are the dim names consistent with cice history output?
    assert set(grids["cice"].ds.dims) == set(
        ["ni", "nj"]
    ), "cice dimension names should be 'ni','nj' to be consistent with history output"
    assert grids["cice"].ds.sizes["ni"] == len(grids["test_ds"].nx)
    assert grids["cice"].ds.sizes["nj"] == len(grids["test_ds"].ny)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_cice_grid(grids):
    # Test : Is the data the same as the test_grid
    for jVar in grids["test_ds"].variables:
        assert_allclose(
            grids["cice"].ds[jVar], grids["test_ds"][jVar], rtol=1e-13, verbose=True, err_msg=f"{jVar} mismatch"
        )


def test_cice_kmt(mom_grid, grids):
    # Test : does the mask match
    mask = mom_grid.mask_ds.mask
    kmt = grids["cice"].kmt_ds.kmt

    assert_allclose(mask, kmt, rtol=1e-13, verbose=True, err_msg="mask mismatch")


def test_cice_grid_attributes(grids):
    # Test: do the expected attributes to exist in the cice ds
    # To-do: rewrite test using the CF-checker (or similar)
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
        "anglet": {
            "standard_name": "angle_of_rotation_from_east_to_x",
            "units": "radians",
            "grid_mapping": "crs",
            "coordinates": "tlat tlon",
        },
        "htn": {"units": "cm", "coordinates": "ulat tlon", "grid_mapping": "crs"},
        "hte": {"units": "cm", "coordinates": "tlat ulon", "grid_mapping": "crs"},
    }

    for iVar in grids["cice"].ds.keys():
        if iVar != "crs":  # test seperately
            for jAttr in cf_attributes[iVar].keys():
                assert grids["cice"].ds[iVar].attrs[jAttr] == cf_attributes[iVar][jAttr]


def test_crs_exist(grids):
    # Test: has the crs been added ?
    # todo: open with GDAL and rioxarray and confirm they find the crs?
    assert hasattr(grids["cice"].ds, "crs")
    assert hasattr(grids["cice"].kmt_ds, "crs")


def test_inputs_logged(grids, mom_grid):
    # Test: have the source data been logged ?

    input_md5 = run(["md5sum", mom_grid.path], capture_output=True, text=True)
    input_md5 = input_md5.stdout.split(" ")[0]
    mask_md5 = run(["md5sum", mom_grid.mask_path], capture_output=True, text=True)
    mask_md5 = mask_md5.stdout.split(" ")[0]

    for ds in [grids["cice"].ds, grids["cice"].kmt_ds]:
        assert ds.inputfile == (
            mom_grid.path + " (md5 hash: " + input_md5 + "), " + mom_grid.mask_path + " (md5 hash: " + mask_md5 + ")"
        ), "inputfile attribute incorrect ({ds.inputfile} != {mom_grid.path})"

        assert hasattr(ds, "inputfile"), "inputfile attribute missing"

        assert hasattr(ds, "history"), "history attribute missing"


def test_variant(mom_grid, tmp_path):
    # Is a error given for variant not equal to None or 'cice5-auscom'

    mom = MomGrid.fromfile(mom_grid.path, mask_file=mom_grid.mask_path)
    cice = CiceGrid.fromgrid(mom)

    # invalid variant (="andrew")
    with pytest.raises(NotImplementedError, match="andrew not recognised"):
        cice.write(str(tmp_path) + "/grid2.nc", str(tmp_path) + "/kmt2.nc", variant="andrew")

    # valid variant (="cice5-auscom")
    try:
        cice.write(str(tmp_path) + "/grid2.nc", str(tmp_path) + "/kmt2.nc", variant="cice5-auscom")
    except:
        assert False, "Failed to write cice grid with valid input arguments provided"

    # valid variant (default = None)
    try:
        cice.write(str(tmp_path) + "/grid2.nc", str(tmp_path) + "/kmt2.nc")
    except:
        assert False, "Failed to write cice grid with 'None' variant"
