#
# R Script: PFM to JASPAR Format Converter
# =========================================
#
# DESCRIPTION:
# This script converts a Position Frequency Matrix (PFM) from a tab-separated file
# to the JASPAR format. The JASPAR format.
#
# RELATED TO compute_motif
# 
# USAGE:
# ```bash
# Rscript pfm_to_jaspar.R <input_matrix> <output_file>
# ```
#
# ARGUMENTS:
#  `<input_matrix>`  : Path to the input PFM file. Required.           
#  `<output_file>`   : Path to the output JASPAR format file. Required.
#
# INPUT FORMAT:
# The input file should be a tab-separated text file where:
# - Each row represents a position in the motif.
# - Each column represents counts for a specific nucleotide.
# - The first row may contain metadata or headers (which will be removed).
#
# Example Input:
# A    C    G    T
# 2    0    1    3
# 1    3    2    0
# 0    2    1    3
#
# OUTPUT FORMAT:
# The output file will be in JASPAR format where:
# - Each row represents a position in the motif.
# - Columns are labeled A, C, G, T.
# - Values represent nucleotide counts at each position.
#
# Example Output:
# A    C    G    T
# 1    3    2    0
# 0    2    1    3
#
# NOTES:
# - The script assumes the input file is properly formatted as a PFM.
# - The first row of the input is skipped, assuming it contains headers or metadata.
# - The output file will have the same number of rows as the input minus one.
#

# LOAD ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
matrix <- args[1]
output<- args[2]

# LOAD PFM MATRIX
pfm <- read.table(matrix,header=FALSE, sep="\t")

# TRANSPOSE
jaspar_pfm <- t(pfm)
colnames(jaspar_pfm) <- c("A","C","G","T")
jaspar_pfm <- jaspar_pfm[2:dim(jaspar_pfm)[1],]
write.table(jaspar_pfm, file=output,quote=FALSE,sep="\t",row.names=FALSE)
