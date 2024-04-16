from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("esmgrids")
except PackageNotFoundError:
    # package is not installed
    pass
