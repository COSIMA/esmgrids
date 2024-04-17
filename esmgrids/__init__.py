import re
import warnings
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("esmgrids")
except PackageNotFoundError:
    # package is not installed
    pass


def safe_version():
    """
    Returns the version, issuing a warning if there are revisions since the last tag
    and an error if there are uncommitted changes

    This function assumes the setuptools_scm default versioning scheme - see
    https://setuptools-scm.readthedocs.io/en/latest/usage/#default-versioning-scheme
    """
    if re.match(r".*\d{8}$", __version__):
        warnings.warn(
            (
                "There are uncommitted changes! Commit, push and release these changes before "
                "generating any production files."
            )
        )
    elif re.match(r".*dev.*", __version__):
        warnings.warn(("There are unreleased commits! Do a release before generating any production files."))

    return __version__