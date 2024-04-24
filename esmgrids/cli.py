import os
import argparse

from esmgrids import safe_version
from esmgrids.util import md5sum
from esmgrids.mom_grid import MomGrid
from esmgrids.cice_grid import CiceGrid


def cice_from_mom():
    """Script for creating CICE grid files from MOM grid files"""

    parser = argparse.ArgumentParser(description="Create CICE grid files from MOM grid files")
    parser.add_argument("--ocean_hgrid", type=str, help="Input MOM ocean_hgrid.nc supergrid file")
    parser.add_argument("--ocean_mask", type=str, help="Input MOM ocean_mask.nc mask file")
    parser.add_argument("--cice_grid", type=str, default="grid.nc", help="Output CICE grid file")
    parser.add_argument("--cice_kmt", type=str, default="kmt.nc", help="Output CICE kmt file")
    parser.add_argument("--cice_variant", type=str, default=None, help="Cice variant")

    args = parser.parse_args()
    ocean_hgrid = os.path.abspath(args.ocean_hgrid)
    ocean_mask = os.path.abspath(args.ocean_mask)
    cice_grid = os.path.abspath(args.cice_grid)
    cice_kmt = os.path.abspath(args.cice_kmt)
    cice_variant = args.cice_variant

    version = safe_version()
    runcmd = (
        f"Created using https://github.com/COSIMA/esmgrids {version}: "
        f"cice_from_mom --ocean_hgrid={ocean_hgrid} --ocean_mask={ocean_mask} "
        f"--cice_grid={cice_grid} --cice_kmt={cice_kmt} --cice_variant={cice_variant}"
    )
    provenance_metadata = {
        "inputfile": (
            f"{ocean_hgrid} (md5 hash: {md5sum(ocean_hgrid)}), " f"{ocean_mask} (md5 hash: {md5sum(ocean_mask)})"
        ),
        "history": runcmd,
    }

    mom = MomGrid.fromfile(ocean_hgrid, mask_file=ocean_mask)
    cice = CiceGrid.fromgrid(mom)
    cice.write(cice_grid, cice_kmt, metadata=provenance_metadata, variant=cice_variant)
