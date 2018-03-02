#!/bin/bash
# Execute in one of the phases directory with the following
# command: bash evol_ft.sh number_of_steps
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Name of the output file with the forces
outfile="evol_ft.csv"

# Initializating and writing header to output file
echo "timestep,total resulting force,horizontal plane resulting force" > ${outfile}

# Number of steps to be passed as an argument
b=$1
for (( c=1; c<=$b; c++ ))
do
	## Computing the total force on all the focal adhesions at a given timestep
	python ${SCRIPTS_DIR}/pfocals_results.py $c ${outfile}
done
