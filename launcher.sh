#!/bin/bash
# Script generating the configuration files for the two stages of the adhesion
# simulation. The two stages consist in a spreading phase of the cell (SPRD) and
# a stabilisation phase of the cell (STBL)

# Time-steping and output frequency parameters of the two phases
## Spreading
nsteps_sprd=1000
dt_sprd=0.0000005
noutp_sprd=20
ndisp_sprd=5
## Stabilization
nsteps_stbl=4000
dt_stbl=0.00000025
noutp_stbl=20
ndisp_stbl=20

# Generation of the configuration files for the two phases of the simulation
python ./gen_sample.py $nsteps_sprd $dt_sprd $nsteps_stbl $dt_stbl $noutp_sprd

# Simulation of the spreading phase
echo "###  PHASE 1 - SPREADING OF THE CELL ###"
cp -r DATBOX_SPRD DATBOX
python ./command.py $nsteps_sprd $dt_sprd $noutp $ndisp_sprd
rm -r phase_1
mkdir phase_1/
cp -r DATBOX phase_1/
cp -r DISPLAY phase_1/
cp -r OUTBOX phase_1/
cp -r POSTPRO phase_1/
rm -r DISPLAY POSTPRO OUTBOX DATBOX

# Transfering data from the spreading phase results to the stabilization phases
# inputs
cp ./DATBOX_SPRD/BODIES.DAT ./DATBOX_STBL/BODIES.OLD
cp ./phase_1/OUTBOX/Vloc_Rloc.LAST ./DATBOX_STBL/
cp ./phase_1/OUTBOX/DOF.LAST ./DATBOX_STBL/
cd DATBOX_STBL
python ../scripts/bodies_dof_update.py
python ../scripts/clean_vloc_dof.py
cd ../

# Simulation of the stabilisation phase
echo "###  PHASE 2 - STABILIZATION OF THE CELL ###"
cp -r DATBOX_STBL DATBOX
python ./command.py $nsteps_stbl $dt_stbl $noutp_stbl $ndisp_stbl
rm -r phase_2
mkdir phase_2/
cp -r DATBOX phase_2/
cp -r DISPLAY phase_2/
cp -r OUTBOX phase_2/
cp -r POSTPRO phase_2/
rm -r DISPLAY POSTPRO OUTBOX DATBOX
