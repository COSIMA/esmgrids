
from __future__ import print_function

import sys
import os
import numpy as np
import numba
import tempfile
import subprocess as sp
import netCDF4 as nc
from scipy import ndimage as nd


@numba.jit
def apply_weights(src, dest_shape, n_s, n_b, row, col, s):
    """
    Apply ESMF regirdding weights.
    """

    dest = np.ndarray(dest_shape).flatten()
    dest[:] = 0.0
    src = src.flatten()

    for i in xrange(n_s):
        dest[row[i]-1] = dest[row[i]-1] + s[i]*src[col[i]-1]

    return dest.reshape(dest_shape)


def regrid(regrid_weights, src_data, dest_grid):
    """
    Regrid a single time index of data.
    """

    print('Horizontal regridding ...')
    # Destination arrays
    dest_data = np.ndarray((dest_grid.num_lat_points,
                            dest_grid.num_lon_points))

    with nc.Dataset(regrid_weights) as wf:
        n_s = wf.dimensions['n_s'].size
        n_b = wf.dimensions['n_b'].size
        row = wf.variables['row'][:]
        col = wf.variables['col'][:]
        s = wf.variables['S'][:]

    dest_data[:, :] = apply_weights(src_data[:, :], dest_data.shape,
                                    n_s, n_b, row, col, s)
    return dest_data


def create_regrid_weights(src_grid, dest_grid, method,
                          unmasked_src=True, unmasked_dest=False):

    assert method == 'bilinear' or method == 'neareststod'

    _, src_grid_scrip = tempfile.mkstemp(suffix='.nc')
    if unmasked_src:
        src_grid.write_scrip(src_grid_scrip,
                             mask=np.zeros_like(src_grid.mask_t, dtype=int))
    else:
        src_grid.write_scrip(src_grid_scrip)

    _, dest_grid_scrip = tempfile.mkstemp(suffix='.nc')
    if unmasked_dest:
        dest_grid.write_scrip(dest_grid_scrip,
                              mask=np.zeros_like(dest_grid.mask_t, dtype=int))
    else:
        dest_grid.write_scrip(dest_grid_scrip)

    _, regrid_weights = tempfile.mkstemp(suffix='.nc')

    try:
        sp.check_output(['ESMF_RegridWeightGen', '-s', src_grid_scrip,
                         '-d', dest_grid_scrip, '-m', method,
                         '-w', regrid_weights])
    except sp.CalledProcessError as e:
        fstring = "Error: ESMF_RegridWeightGen failed return code {}"
        print(fstring.format(e.returncode), file=sys.stderr)
        print(e.output, file=sys.stderr)
        log = 'PET0.RegridWeightGen.Log'
        if os.path.exists(log):
            print('Contents of {}:'.format(log), file=sys.stderr)
            with open(log) as f:
                print(f.read(), file=sys.stderr)
        return None

    assert(os.path.exists(regrid_weights))

    return regrid_weights
