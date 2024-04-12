import numpy as np
import pyproj
from shapely.geometry import shape
import subprocess
import os
from datetime import datetime
from netCDF4 import Dataset


proj_str = "+proj=laea +lat_0={} +lon_0={} +ellps=sphere"


def calc_area_of_polygons(clons, clats):
    """
    Calculate the area of lat-lon polygons.

    We project sphere onto a flat surface using an equal area projection
    and then calculate the area of flat polygon.

    This is slow we should do some caching to avoid recomputing.

    FIXME: compare against a Haversine method for speed and accuracy
    """

    areas = np.zeros(clons.shape[1:])
    areas[:] = np.NAN

    for j in range(areas.shape[0]):
        for i in range(areas.shape[1]):

            lats = clats[:, j, i]
            lons = clons[:, j, i]

            lat_centre = lats[0] + abs(lats[2] - lats[1]) / 2
            lon_centre = lons[0] + abs(lons[1] - lons[0]) / 2

            pa = pyproj.Proj(proj_str.format(lat_centre, lon_centre))
            x, y = pa(lons, lats)

            cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
            areas[j, i] = shape(cop).area

    assert np.sum(areas) is not np.NAN
    assert np.min(areas) > 0

    return areas


def is_git_repo():
    """
    Return True/False depending on whether or not the current directory is a git repo.
    """

    return subprocess.call(["git", "-C", ".", "status"], stderr=subprocess.STDOUT, stdout=open(os.devnull, "w")) == 0


def git_info():
    """
    Return the git repo origin url, relative path to this file, and latest commit hash.
    """

    url = subprocess.check_output(["git", "remote", "get-url", "origin"]).decode("ascii").strip()
    top_level_dir = subprocess.check_output(["git", "rev-parse", "--show-toplevel"]).decode("ascii").strip()
    rel_path = os.path.relpath(__file__, top_level_dir)
    hash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()

    return url, rel_path, hash


def create_nc(filename):

    f = Dataset(filename, "w")

    f.timeGenerated = f"{datetime.now()}"
    f.created_by = f"{os.environ.get('USER')}"
    if is_git_repo():
        git_url, _, git_hash = git_info()
        f.history = f"Created using commit {git_hash} of {git_url}"

    return f


def md5sum(filename):
    from hashlib import md5
    from mmap import mmap, ACCESS_READ

    with open(filename) as file, mmap(file.fileno(), 0, access=ACCESS_READ) as file:
        return md5(file).hexdigest()
