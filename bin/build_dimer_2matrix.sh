#!/bin/bash
# This script generates a dimer matrix by combining two input matrices based on the specified dimer type and spacing.
# 
# RELATED TO: spacing_2TFs
# 
# USAGE:
#   ./build_dimer_2matrix.sh <matrix_1> <offset_left1> <offset_right1> <matrix_2> <offset_left2> <offset_right2> <name_matrix> <dimer_type> <spacing>
#
# ARGUMENTS:
#   matrix_1      - Path to the first input matrix file.
#   offset_left1  - Left offset for the first matrix.
#   offset_right1 - Right offset for the first matrix.
#   matrix_2      - Path to the second input matrix file.
#   offset_left2  - Left offset for the second matrix.
#   offset_right2 - Right offset for the second matrix.
#   name_matrix   - Name for the resulting dimer matrix.
#   dimer_type    - Type of dimer to generate. Supported types:
#                   - IR: Inverted Repeat (-> <-)
#                   - DR: Direct Repeat (-> ->)
#                   - ER: Everted Repeat (<-->)
#                   - DR12/DRab: Direct Repeat with specific offsets.
#                   - DR21/DRba: Direct Repeat with reversed offsets.
#   spacing       - Spacing between the two matrices in the dimer.
#
# OUTPUT:
#   The script outputs the resulting dimer matrix to the standard output.
#
# FUNCTIONALITY:
#   - Reads the input matrices and their respective offsets.
#   - Calculates the spacing and offsets based on the dimer type.
#   - Combines the matrices according to the specified dimer type:
#     - IR: Inverts the second matrix and appends it after the first matrix with spacing.
#     - DR: Appends the second matrix after the first matrix with spacing.
#     - ER: Inverts the first matrix and appends it before the second matrix with spacing.
#     - DR12/DRab: Similar to DR, with specific offset handling.
#     - DR21/DRba: Similar to DR, but the order of matrices is reversed.
#   - Outputs the resulting dimer matrix in the specified format.

# FUNCTION
calc(){
awk "BEGIN { print "$*" }";
}

# LOAD ARGUMENTS
matrix_1=$1
offset_left1=$2
offset_right1=$3
matrix_2=$4
offset_left2=$5
offset_right2=$6
name_matrix=$7
dimer_type=$8
spacing=$9

header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE\nA\tC\tG\tT"

pfm1_mono=$(tail -n +3 $matrix_1)
pfm2_mono=$(tail -n +3 $matrix_2)

if [ $dimer_type == "IR" ] #-> <- IR
    then
    offset=$(calc $offset_right1+$offset_right2)
    spacing=$(calc $spacing-$offset)
    echo -e "$header"
    echo -e "$pfm1_mono"
    for i in `seq 1 $spacing`; do echo -e "1\t1\t1\t1"; done
    echo -e "$pfm2_mono" | tac | awk -v OFS="\t" '{print $4,$3,$2,$1}'
elif [ $dimer_type == "DR" ] # ->-> DR
    then
    offset=$(calc $offset_right1+$offset_left2)
    spacing=$(calc $spacing-$offset)
    echo -e "$header"
    echo -e "$pfm1_mono"
    for i in `seq 1 $spacing`; do echo -e "1\t1\t1\t1"; done
    echo -e "$pfm2_mono"
elif [ $dimer_type == "ER" ] # ER <-->
    then
    offset=$(calc $offset_left1+$offset_left2)
    spacing=$(calc $spacing-$offset)
    echo -e "$header"
    echo -e "$pfm1_mono" | tac | awk -v OFS="\t" '{print $4,$3,$2,$1}'
    for i in `seq 1 $spacing`; do echo -e "1\t1\t1\t1"; done
    echo -e "$pfm2_mono"
elif  [ $dimer_type == "DR12" ] || [ $dimer_type == "DRab" ] # ->-> DR
    then
    offset=$(calc $offset_right1+$offset_left2)
    spacing=$(calc $spacing-$offset)
    echo -e "$header"
    echo -e "$pfm1_mono"
    for i in `seq 1 $spacing`; do echo -e "1\t1\t1\t1"; done
    echo -e "$pfm2_mono"
elif  [ $dimer_type == "DR21" ] || [ $dimer_type == "DRba" ] # ->-> DR
    then
    offset=$(calc $offset_right2+$offset_left1)
    spacing=$(calc $spacing-$offset)
    echo -e "$header"
    echo -e "$pfm2_mono"
    for i in `seq 1 $spacing`; do echo -e "1\t1\t1\t1"; done
    echo -e "$pfm1_mono"
else
    echo "incorrect dimer, dimer type specified (ER,IR,DR)" >&2
    exit 1
fi
exit 0

# ##############################################################################################################
#### AFTER THIS POINT: OLD version, supporting only same matrix

matrix=$1
name_matrix=$2
spacing=$3
dimer_type=$4
offset_left=$5
offset_right=$6

calc(){
awk "BEGIN { print "$*" }";
}


if [ -z $dimer_type ]
then
    echo "dimer type is  ER, IR, DR">&2
    exit 1
fi


header="MATRIX COUNT SYMMETRIC $name_matrix SIMPLE\nA\tC\tG\tT"

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
    header="MATRIX COUNT SIMPLE $name_matrix SIMPLE\nA\tC\tG\tT"
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


