#
# R Script: Sequence Logo Generator
# =================================
#
# DESCRIPTION:
# This script generates a sequence logo from SeqConv generated matrices.
# 
# RELATED TO compute_SeqConv
# 
# USAGE:
# ```bash
# Rscript sequence_logo.R <pfm_file> <output_file>
# ```
#
# ARGUMENTS
# - `<pfm_file>`       : Path to the input PFM file. Required.           
# - `<output_file>`    : Path to save the output PNG image. Required.   
# 
# INPUT FORMAT:
# The input PFM file should be a tab-separated text file where:
# - Each row represents a position in the motif.
# - Each column represents the frequency of a nucleotide (A, C, G, T).
# - The file should contain exactly 4 columns (one for each nucleotide).
#
# Example Input:
# 0.2 0.3 0.1 0.4
# 0.1 0.4 0.2 0.3
# 0.3 0.1 0.4 0.2
#
# OUTPUTS:
# The logo visually represents the frequency of each nucleotide at each position
# in the motif, with letter height proportional to frequency and color-coded by nucleotide.
# 
# NOTES:
# - The script assumes the input PFM contains valid frequency values (0-1).
# - The script automatically adjusts the x-axis limit to 24 positions, which
#   should be updated if working with motifs of different lengths.
#

# LOAD LIBRARIES
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/", .libPaths()))
library(ggseqlogo)
library(ggplot2)

# ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)

# LOAD DATA 
pfm_data <- read.table(args[1], sep = "\t", header = FALSE)

# Transpose matrix
pfm_matrix <- t(as.matrix(pfm_data))
# Add appropriate row and column names
rownames(pfm_matrix) <- c("A", "C", "G", "T")
colnames(pfm_matrix) <- paste0("Pos", 1:ncol(pfm_matrix))
# Create a sequence logo object
colors <- make_col_scheme(
  chars = c('A', 'C', 'G', 'T'),
  cols = c('red', 'blue', 'yellow', 'green'),
  name = 'colors'
)
seq_data <- ggseqlogo(pfm_matrix, seq_type = "dna", col_scheme = colors)
seq_data <- seq_data + expand_limits(x = 24, y = 2)

# SAVE AND EXPORT 
# print
print(seq_data)
# Save the plot as a PNG file
ggsave(args[2], width = 10, height = 5, dpi = 300)

