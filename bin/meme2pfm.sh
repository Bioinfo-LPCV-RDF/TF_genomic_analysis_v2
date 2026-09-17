#!/bin/bash

#
# Script: MEME Matrix Converter
# =============================
#
# DESCRIPTION
# This script converts MEME motif files into separate PFM (Position Frequency Matrix) files.
# It processes a MEME output file, extracts individual motifs, and converts them into PFM format.
# The script handles both asymmetric and symmetric motifs based on the provided type.
#
# RELATED TO compute_motif
# 
# USAGE:
# ```bash
# ./script.sh <meme_matrix_file> <name_matrix> <type_pal>
# ```
#
# ARGUMENTS:
# - `<meme_matrix_file>`  : Path to the MEME output file containing motifs. Required.               
# - `<name_matrix>`       : Base name for the output matrices. Required.                                
# - `<type_pal>`          : Type of motif: "asym" for asymmetric or "sym" for symmetric. Required.     
#
# NOTES:
# ------
# - The script assumes the MEME file follows the standard MEME format.
# - The output PFM files are saved in a subdirectory of the input file's directory.
# - The script handles both symmetric and asymmetric motifs based on the type_pal parameter.
# - The PFM format includes a header line followed by the matrix of nucleotide counts.
#
# OUTPUTS
# Each PFM file contains:
# 1. A header line indicating the matrix type (symmetric/asymmetric) and name.
# 2. A matrix of nucleotide counts (A, C, G, T) for each position in the motif.
#
# Example PFM File:
# MATRIX COUNT ASYMMETRIC Motif_1 SIMPLE
# A   C   G   T
# 10  2   3   5
# 5   8   7   0
# ...
#

# LOAD ARGUMENTS
meme_matrix=$1
name_matrix=$2
type_pal=$3

# PATH HANDLING
output_path=$(echo $meme_matrix |awk '{ gsub("/[/]*","/",$1);split($1,PATH,"/"); if(length(PATH) >1) {if(PATH[1] != ""){ OUT_PATH=PATH[1]; for (i=2; i<length(PATH); i++){OUT_PATH=OUT_PATH"/"PATH[i]} } else {OUT_PATH="/"; for (i=2; i<length(PATH); i++){OUT_PATH=OUT_PATH"/"PATH[i]}  } } else {OUT_PATH="./"} } END{print OUT_PATH}')
# get start line
start_line=$(awk 'BEGIN{FS="\t"; OCCUR=1} $0 ~ "MOTIF" && OCCUR!=1{VAR=VAR"\t"FNR} $0 ~ "MOTIF" && OCCUR==1{VAR=FNR; OCCUR=2} END{VAR=VAR"\t"FNR+2; print VAR}' $meme_matrix)
mkdir -p -m 774 ${output_path}/Motif_MEME_seperateFiles

# PFM CONVERSION
echo $start_line | awk 'BEGIN{MOTIF_NUMBER=1} NR==FNR{for (i=1; i<=NF; i++){TABLE[i]=$i} } NR!=FNR{if(FNR>=TABLE[MOTIF_NUMBER] && FNR<TABLE[MOTIF_NUMBER+1]-1){print $0 > "'${output_path}'/Motif_MEME_seperateFiles/Motif_"MOTIF_NUMBER".txt"} else if(FNR == TABLE[MOTIF_NUMBER+1]-1){print $0 "'${output_path}'/Motif_MEME_seperateFiles/Motif_"MOTIF_NUMBER".txt" ; MOTIF_NUMBER+=1} } ' - $meme_matrix

for motif_file in $(ls ${output_path}/Motif_MEME_seperateFiles/Motif_*.txt)
do
	name_matrix=$(head -1 $motif_file | awk '{print $3"_"$2}')
	name_file=$(head -1 $motif_file | awk '{print $3}' | awk -v FS="-" '{print "Motif_"$2}')
	if [ $type_pal == "asym" ]; then
		header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE"
	else
		header="MATRIX COUNT SYMMETRIC $name_matrix SIMPLE"
	fi
# 	header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE"
	matrix=$(tail -n +3 $motif_file)
	line1=$(head -1 <(echo -e "$matrix"))
	line1=( $line1 )
	nsites=${line1[7]}
	pfm=$(tail -n +2 <(echo -e "$matrix") | sed -e 's/^\s\s*//' -e 's/\s\s*/\t/g' -e 's/\s\s*$//g' | awk -v n=$nsites -v OFS="\t" '{print int($1*n)+1,int($2*n)+1,int($3*n)+1,int($4*n)+1}'| sed "1iA\tC\tG\tT")
	echo -e "$header\n$pfm" > $(dirname ${motif_file})"/"${name_file}".pfm"
	rm ${motif_file}
done

# TODO: to rm ? 
# #!/bin/bash
# 
# meme_matrix=$1
# name_matrix=$2
# 
# header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE"
# 
# start_line=$(grep -n "letter-probability matrix" $meme_matrix | awk -v FS=":" '{print $1}')
# matrix=$(tail -n +$start_line $meme_matrix)
# line1=$(head -1 <(echo -e "$matrix"))
# line1=( $line1 )
# nsites=${line1[7]}
# pfm=$(tail -n +2 <(echo -e "$matrix") | sed -e 's/^\s\s*//' -e 's/\s\s*/\t/g' -e 's/\s\s*$//g' | awk -v n=$nsites -v OFS="\t" '{print int($1*n)+1,int($2*n)+1,int($3*n)+1,int($4*n)+1}'| sed "1iA\tC\tG\tT")
# echo -e "$header\n$pfm"
# 
# exit 0