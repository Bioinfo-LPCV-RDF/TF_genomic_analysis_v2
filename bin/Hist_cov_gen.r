# RScript: Coverage Ratio Comparison

## DESCRIPTION 
# This R script generates scatter plots to compare coverage ratios between two samples. 
# It visualizes the relationship between the coverage of two datasets and highlights the coverage fold change (CFC). 
# The script uses ggplot2 and cowplot to create high-quality, publication-ready plots.

## RELATED TO pairwize_comparison

## USAGE
# The script is designed to be called from the command line with the following arguments:
# ```bash
# Rscript Hist_cov_gen.r <output_directory> <sample_name1> <sample_name2> <table_file> <type>
## ARGUMENTS
# - `<output_directory>`: The directory where the output plots will be saved.
# - `<sample_name1>`: The name of the first sample to compare.
# - `<sample_name2>`: The name of the second sample to compare.
# - `<table_file>`: The path to the input table file (tab-delimited) containing coverage data.
# - `<type>`: The type of comparison (not used in the current script).

## OUTPUT
# <sample1_name>_<sample2_name>.png: A combined plot showing the relationship between the coverage 
# of the two samples, with legends for peak calling and coverage fold change.

## NOTES
# - Logarithmic Scales: Both axes use a logarithmic scale to better visualize the data distribution.

# LOAD LIBRARIES
rm(list=ls())
# library(Biostrings ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(ggplot2 ,verbose=FALSE ,warn.conflicts=FALSEd, quietly = TRUE,lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(Cairo ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE,lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/") 
library(cowplot,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE,lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")

# LOAD ARGUMENTS
args=commandArgs(trailingOnly=TRUE)
out_dir=args[1]
name1=args[2]
name2=args[3]
tableName=args[4]
Type=args[5]

# LOAD DATA
# Rscript Hist_cov_gen.r dir name1 name2 tablename rationame1 rationame2
dir.create(file.path(out_dir),recursive=TRUE ,showWarnings = FALSE)
table <- read.table(tableName,sep="\t",header=TRUE,check.names=FALSE)

table$ratio<-table[,name1] / table[,name2]
table$name <- factor(table$name, levels = c(name1,name2,"both"))

head(table)

# PLOTTING
# MACS3 peaks attribution plot
covRiL1bis<-ggplot(data=table,aes(y=table[,name1],x=table[,name2]))+
geom_point(alpha=0.4,aes(color=name)) +
scale_x_log10() +
scale_y_log10() +
labs(color="peaks called \nby MACS3") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))

# CFC peaks attribution plot
covRiL2bis<-ggplot(data=table,aes(y=table[,name1],x=table[,name2]))+
geom_point(alpha=0.4,aes(color=log10(ratio))) +
scale_x_log10() +
scale_y_log10() +
scale_color_gradientn(colours=rainbow(6),name=paste("CFC = \n", name1," / ",name2,sep=""))+
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))

legend1bis <-get_legend(covRiL1bis)
covRiL1bis <- covRiL1bis + theme(legend.position='none')

legend2bis <-get_legend(covRiL2bis)
covRiL2bis <- covRiL2bis + theme(legend.position='none')

# export
Cairo(width = 1500, height = 500, file=paste(out_dir,"/",name1,"_",name2,".png",sep=""), type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
# plot_grid(covRiL1, covRiL2, nrow=1, ncol=2)
ggdraw() +
draw_plot(covRiL1bis, 0.0, 0.0, 0.5, 1) +
draw_plot(covRiL2bis, 0.5, 0, 0.5, 1) +

draw_plot(legend1bis, 0.43, -0.02, 0.2, 0.5)+
draw_plot(legend2bis, 0.90, -0.02, 0.18, 0.5)+

draw_plot_label(c("Reads mapped in library", "Reads mapped in library"), c(0.02, 0.52), c(0.995, 0.995), size = 15)

dev.off()

# }


