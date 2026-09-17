#
# R Script: Score Distribution Analyzer
# ====================================
#
# DESCRIPTION
# This R script analyzes the distribution of scores from a dataset, calculates
# percentiles, and generates visualizations.
#
# The script reads a tab-separated file containing scores, calculates various
# percentiles (especially in the 99th percentile range), generates a density plot,
# and saves the results to files.
# 
# RELATED TO compute_distribution
#
# USAGE
# The script is designed to be called from the command line with the following arguments:
# ```bash
# Rscript score_analyzer.R <input_file> <output_dir> <percentile>
# ```
#
# ARGUMENTS
#  `<input_file>`    : Path to the input tab-separated file containing scores. Required. 
#  `<output_dir>`    : Directory where output files will be saved. Required. 
#  `<percentile>`    : Percentile value to calculate (e.g., 0.995). Required. 
#
# INPUT FORMAT 
# The input file should be a tab-separated text file where:
# - Column 8 contains the scores to analyze
# - The file should have a header row (though header=FALSE is used in read.table)
#
# Example Input:
# seq1   pos1   ...   0.85
# seq2   pos2   ...   0.92
# ...
#
# OUTPUTS
# The script generates the following files in the specified output directory:
# - distribution.png: Density plot of score distribution
# - scores_info.txt: Table of percentile values from 99.75% to 99.99%
# - scores_99_0_info.txt: 99.0th percentile value
# - scores_99_20_info.txt: 99.2th percentile value
# - scores_99_50_info.txt: 99.5th percentile value
# - scores_99_70_info.txt: 99.7th percentile value
# - scores_99_90_info.txt: 99.9th percentile value
#
# NOTES
# - The script focuses on high percentiles (99th percentile and above), which is
#   useful for identifying significant outliers or extreme values in the data.
# - The ECDF calculation shows what percentage of scores fall below the specified percentile.
#

# LOAD LIBRARIES
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(ggplot2)
library(dplyr)
library("cowplot",verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)

# LOAD ARGS
args <- commandArgs(TRUE)

# LOAD DATA 
print("[INFO] - Loading dataset")
dataset<-read.table(args[1],sep="\t",header=FALSE)

# PLOT
p<-ggplot(dataset,aes(x=V8))+
#   geom_histogram(aes(y=..density..), fill="lightblue", bins = 60)+
  geom_density(aes(y=..count../sum(..count..)))+
#   stat_ecdf(geom="point")+
  xlab("score")+
#   scale_y_continuous(name = "%", labels=scales::percent)+
  theme_classic()
  save_plot(paste(args[2],"/distribution.png",sep=""), p, base_width=8, base_height=4)
# p <- ggplotly(p)



scores<-dataset[,8]
print("[INFO] - Table read")
## calculate 99th and 95th percentile of scores
perc990_value<-quantile(scores, .990)
perc992_value<-quantile(scores, .992)
perc995_value<-quantile(scores, .995)
perc997_value<-quantile(scores, .997)
perc999_value<-quantile(scores, .999)

# LOOP IN SCORES 
perc_values<-c()
colnam<-c()
i=750
while(i<1000){
perc_values<-c(perc_values,quantile(scores, .9+(i/10000)))
colnam<-c(colnam,.9+(i/10000))
	i=i+1
}

# for each percentage (and score associated), used the length of genome tested to estimate number of scores (a.k.a TFBS) higher than thresholds scores to obtain TFBS/kb
genome_length<- as.numeric(args[3])
# col4 => strand, we keep only scores on + strand for this analysis
scores<-dataset[dataset[,4]=="+",8]

tfbs_values<-c()
for(i in 1:length(perc_values)){
	# check number of scores higher than threshold
	n_scores_higher<-length(scores[scores>=perc_values[i]])
	# estimate number of TFBS per kb
	TFBS_per_kb<-n_scores_higher/(genome_length/1000)
	tfbs_values<-c(tfbs_values,TFBS_per_kb)
}

table<-data.frame(perc_values,colnam,tfbs_values)
write.table(table,file=paste0(args[2],"/scores_info.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)

# same as previous loop but using nb of TFBS per kb as reference instead of percentage of scores, to obtain the percentage of scores associated with a given TFBS/kb and the threshold score associated with this percentage of scores
# loop TFBS/kb from 0.01 to 5
tfbs_values<-c()
colnam<-c()
perc_values<-c()
i=500
while(i>=1){
	
	# estimate number of TFBS per kb
	TFBS_per_kb<-c(i/100)
	# estimate number of scores higher than threshold
	n_scores_higher<-TFBS_per_kb*(genome_length/1000)
	# check percentage of scores higher than threshold
	perc_score_higher<-n_scores_higher/length(scores)
	# estimate threshold score associated with this percentage of scores
	threshold_score<-quantile(scores, 1-perc_score_higher)
	# register in vectors
	tfbs_values<-c(tfbs_values,TFBS_per_kb)
	colnam<-c(colnam,threshold_score)
	perc_values<-c(perc_values,1-perc_score_higher)
	i=i-1
}
table<-data.frame(colnam,perc_values,tfbs_values)
write.table(table,file=paste0(args[2],"/scores_info_TFBS_kb.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)






# PLOT AND WRITE OUT TABLE 
ecdf_fun <- function(x,perc) ecdf(x)(perc)
print("[INFO] - Requested percentage for score input is:")

write.table(paste("99_0 percentile:",perc990_value,sep="\t"),
			file=paste0(args[2],"/scores_99_0_info.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)
write.table(paste("99_5 percentile:",perc995_value,sep="\t"),
			file=paste0(args[2],"/scores_99_50_info.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)
write.table(paste("99_9 percentile:",perc999_value,sep="\t"),
			file=paste0(args[2],"/scores_99_90_info.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)
write.table(paste("99_2 percentile:",perc992_value,sep="\t"),
			file=paste0(args[2],"/scores_99_20_info.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)
write.table(paste("99_7 percentile:",perc997_value,sep="\t"),
			file=paste0(args[2],"/scores_99_70_info.txt"),quote=FALSE,
			row.names = FALSE,col.names = FALSE)
