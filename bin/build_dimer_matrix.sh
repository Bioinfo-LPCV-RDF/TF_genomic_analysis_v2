#!/bin/bash

# Dimeric Motif Matrix Generator
# ==============================

# DESCRIPTION
# This Bash script generates Position Frequency Matrices (PFM) for homo-dimeric DNA motifs
# from a single monomer matrix. It creates different types of dimers (IR, DR, ER) based
# on provided parameters.

# RELATED TO: compute_ROCS

# Supported Dimer Types:
# - IR (Inverted Repeat): (-> <-)
# - DR (Direct Repeat): (-> ->)
# - ER (Everted Repeat): (<->)

# ARGUMENTS
# The script takes the following parameters in order:
# 1. matrix       : Path to the monomer matrix file
# 2. name_matrix  : Name for the generated dimeric matrix
# 3. spacing      : Number of base pairs between monomers
# 4. dimer_type   : Type of dimer (ER, IR, DR)
# 5. offset_left  : Left offset for positioning
# 6. offset_right : Right offset for positioning

# INPUT FORMAT 
# The matrix file should be in standard PFM format:
# MATRIX COUNT [SYMMETRIC|ASYMMETRIC] <name> SIMPLE
# A    C    G    T
# <values for position 1>
# ...
# <values for position N>

# NOTES
# - Spacer regions are filled with uniform values (1 for all bases)
# - Offsets allow precise positioning of monomers

# LOAD ARGUMENTS
matrix=$1
name_matrix=$2
spacing=$3
dimer_type=$4
offset_left=$5
offset_right=$6

# FUNCTION
calc(){
awk "BEGIN { print "$*" }";
}

# CHECK
if [ -z $dimer_type ]
then
    echo "[ERROR] - Dimer type is either ER, IR, DR">&2
    exit 1
fi

# INIT
header="MATRIX COUNT SYMMETRIC $name_matrix SIMPLE\nA\tC\tG\tT"

# GENERATE DIMER 
if [ $dimer_type == "IR" ] #-> <- IR
then
    spacing=$(calc $spacing-2*$offset_right)
    pfm_mono=$(tail -n +3 $matrix)
    echo -e "$header"
    echo -e "$pfm_mono"
    for i in `seq 1 $spacing`
    do
	echo -e "1\t1\t1\t1"
    done
    echo -e "$pfm_mono" | tac | awk -v OFS="\t" '{print $4,$3,$2,$1}'
    
elif [ $dimer_type == "DR" ] # ->-> DR
    then
    spacing=$(calc $spacing-$offset_right-$offset_left)
    header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE\nA\tC\tG\tT"
    pfm_mono=$(tail -n +3 $matrix)
    echo -e "$header"
    echo -e "$pfm_mono"
    for i in `seq 1 $spacing`
    do
	echo -e "1\t1\t1\t1"
    done
    echo -e "$pfm_mono" | cat | awk -v OFS="\t" '{print $1,$2,$3,$4}'
    

elif [ $dimer_type == "ER" ] # ER <-->
    then
    spacing=$(calc $spacing-2*$offset_left)
    pfm_mono=$(tail -n +3 $matrix)
    echo -e "$header"
    echo -e "$pfm_mono" | tac | awk -v OFS="\t" '{print $4,$3,$2,$1}'
    for i in `seq 1 $spacing`
    do
	echo -e "1\t1\t1\t1"
    done
    echo -e "$pfm_mono"

echo "$Spacing" >&2
else
    echo "incorrect dimer dimer type specified (ER,IR,DR)" >&2
    exit 1
fi

exit 0
