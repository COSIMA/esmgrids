![Code Health](https://github.com/COSIMA/esmgrids/actions/workflows/ci.yml/badge.svg)

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
pip install '.[tests]'
python -m pytest -m "not broken"
```
