
[![Code Health](https://landscape.io/github/DoublePrecision/esmgrids/master/landscape.svg?style=flat)](https://landscape.io/github/DoublePrecision/esmgrids/master)
[![Build Status](https://travis-ci.org/DoublePrecision/esmgrids.svg?branch=master)](https://travis-ci.org/DoublePrecision/esmgrids)

# esmgrids

This package contains Python representations of Earth System Model grids. These are very useful, for example, for converting model grids between different formats.

Grids currently supported:

- Global regular Lat-Lon
- Global tri-polar
- MOM5 (Modular Ocean Model) tri-polar
- NEMO tri-polar
- GODAS reanalysis
- ORAS4 reanalysis
- CICE tri-polar

## To run the tests

```
conda env create -n grids python3
source activate grids
python -m pytest
```

Warning: this will download a rather large tarball of test inputs.

