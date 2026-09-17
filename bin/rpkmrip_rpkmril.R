#
# R Script: Normalize Read Counts
# ===============================
#
# DESCRIPTION:
# This script normalizes read counts from sequencing experiments.
# It can normalize counts either by the total reads in peaks (inPeaks) or by the total reads in libraries (inLibs).
# The script calculates Counts Per Million (CPM) and Reads Per Kilobase Million (RPKM) for each peak across samples.
#
# RELATED TO compute_rpkmrip_rpkmril, replicates_comaprison, pairwize_comparison
# USAGE:
# ```bash
# Rscript normalize_reads.R <input_file> <label> <working_directory> [sample_names...]
# ```
#
# ARGUMENTS:
#  - `<input_file>`          : Path to the input file containing read counts per peak per sample. A tab-separated file. Required. 
#  - `<label>`               : Normalization mode: "inPeaks" or "inLibs". Required.                        
#  - `<working_directory>`   : Working directory where the script will be executed. Required.          
#  - `[sample_names...]`     : Names of the samples (columns in the input file). Required.                
#
# OUTPUTS:
# - `peaks_perSample_rpkm_[label].txt`: A tab-separated file containing normalized RPKM values.
# 
# NOTES:
# - The script assumes that the input file has at least 4 columns: chromosome, start, end, and read counts for each sample.
# - For "inLibs" mode, the script requires a file named "tmpTotalTags.txt" in the working directory.
# - The output file contains the genomic coordinates (chr, start, end) followed by RPKM values for each sample.
#

# LOAD ARGUMENTS
args <- commandArgs(TRUE)
# print(args)

infile=args[1]
label=args[2] # this is the mode: inPeaks or inLibs
# print(label)

#load the RC per consensus peaks per replicats
d<-read.table(infile, header=F)

# LOAD DATA 
WD=args[3]
setwd(WD)

colnames(d)<-c("chr", "start", "end", args[4:length(args)])
lastCol=ncol(d)

# NORMALIZATION
if (label=="inPeaks"){ 
    # Sum of RC in peaks
    effRC <- colSums(as.matrix(d[,4:lastCol]))  # <- use as.matrix!
    print("-> Total number of reads in peaks per sample:")
    print(effRC)
} else if (label=="inLibs") { 
    # Total number of reads in libraries
    RClibs <- read.table("tmpTotalTags.txt", header = FALSE, stringsAsFactors = FALSE)
    effRC <- suppressWarnings(as.numeric(RClibs$V1))
    effRC <- effRC[!is.na(effRC)]
    print("-> Total number of reads in library per sample ")
    print(effRC)
} else {
    print("[ERROR] - Need to tell the program if normalization should be 'inPeaks' or 'inLibs'")
    quit()
}


#normalize peaks RC by tot RC in consensus filtered peaks or lib size
# 1. Calculates Counts Per Million (CPM) by dividing read counts by the total reads (scaled to 1 million)
CPM<-as.data.frame(t(t(d[,c(4:lastCol)]) / (effRC/1000000)))
print("-> Head of RPM per peak per sample:")
head(CPM)
# 2. Calculates Reads Per Kilobase Million (RPKM) by adjusting CPM for peak length.
RPKM<-as.data.frame(CPM / ((d[,3]-d[,2])/1000))
print("-> Head of RPKM per peak per sample:")
head(RPKM)
print("-> Summary RPKM:")
summary(RPKM)

# FINAL FORMATTING AND EXPORT 
tmp <- cbind(d[,c(1,2,3)],RPKM)
colnames(tmp)<-c("chr", "start", "end", args[4:length(args)])
filename=paste("peaks_perSample_rpkm", label, ".txt", sep="")
write.table(tmp, file=filename, sep="\t", row.names=F, quote=F)

