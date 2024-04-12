#!/bin/bash

#PBS -W umask=0022
#PBS -l mem=24gb
#PBS -l storage=gdata/ik11+gdata/tm70+gdata/hh5
#PBS -l wd
#PBS -j oe

run='python3 esmgrids/cice_grid.py'

umask 0003

module purge
module use /g/data/hh5/public/modules
module load conda/analysis3-23.10

echo "1 degree"
$run /g/data/ik11/inputs/access-om2/input_20201102/mom_1deg/ocean_hgrid.nc /g/data/ik11/inputs/access-om2/input_20201102/mom_1deg/ocean_mask.nc

mkdir 1deg
mv grid.nc kmt.nc 1deg

echo "0.25 deg"
$run /g/data/ik11/inputs/access-om2/input_20230515_025deg_topog/mom_025deg/ocean_hgrid.nc /g/data/ik11/inputs/access-om2/input_20230515_025deg_topog/mom_025deg/ocean_mask.nc

mkdir 025deg
mv grid.nc kmt.nc 025deg

echo "01 deg"
$run /g/data/ik11/inputs/access-om2/input_20201102/mom_01deg/ocean_hgrid.nc /g/data/ik11/inputs/access-om2/input_20201102/mom_01deg/ocean_mask.nc

mkdir 01deg
mv grid.nc kmt.nc 01deg