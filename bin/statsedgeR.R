#' Differential Coverage Analysis with edgeR
#' =========================================
#' 
#' DESCRIPTION
#'  This script performs differential expression analysis using the edgeR package. 
#' It takes as input a count table and outputs a table of results with log fold changes, 
#' p-values, and adjusted p-values (FDR).
#' 
#' USAGE 
#'  Rscript statsedgeR.R <input_table> <output_directory> <group1_size> <group2_size>
#' 
#' ARGUMENTS
#' - input_table Path to the input count table. The table should have columns for 
#'   chromosome, start, end, and counts for each sample.
#' - output_directory Path to the directory where the results will be saved.
#' - group1_size Number of samples in group 1.
#' - group2_size Number of samples in group 2.
#' 
#' PROCESS
#' The script performs the following steps:
#' 1. Reads the input count table and processes it to create unique row names.
#' 2. Defines the experimental groups based on the provided group sizes.
#' 3. Creates a DGEList object and estimates dispersion.
#' 4. Fits a generalized linear model (GLM) and performs a quasi-likelihood F-test.
#' 5. Adjusts p-values using the Benjamini-Hochberg method and calculates fold change categories.
#' 6. Outputs the results to a TSV file in the specified output directory.
#' 
#' OUTPUT
#' A TSV file named `results_edgeR.tsv` containing the following columns:
#' - peaks: Unique identifier for each genomic region.
#' - logFC: Log fold change.
#' - logCPM: Log counts per million.
#' - F: F-statistic.
#' - PValue: Raw p-value.
#' - FDR: Adjusted p-value (False Discovery Rate).
#' - FC: Fold change category (-1, 0, 1).


# LOAD LIBRARIES
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(edgeR)
args=commandArgs(trailingOnly=TRUE)
out_dir=args[2]
tableName=args[1]
nbgrp1=args[3]
nbgrp2=args[4]
print(args)

# rep(2,nbgrp2)
# rep(1,nbgrp1)

data<-read.table(tableName,header = T)
data$peaks<-paste(data$chr,":",data$start,"-",data$end,sep="")
row.names(data)<-make.names(data$peaks, unique=TRUE)
data$chr=NULL
data$start=NULL
data$end=NULL
data$peaks=NULL

group <- factor(c(rep(1,nbgrp1),rep(2,nbgrp2)))
y <- DGEList(counts=data,group=group)
y
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
fit$genes$qlf=decideTests(qlf,adjust.method = "BH", lfc=1)

qlf$table$FDR <- p.adjust(qlf$table$PValue, method="BH")

qlf$table$FC=0
qlf$table$FC[qlf$table$logFC< -1] = -1
qlf$table$FC[qlf$table$logFC> 1] = 1
qlf$table$peaks<-row.names(qlf)
qlf$table<-qlf$table[,c(7,1,2,3,4,5,6)]
write.table(qlf$table,paste(out_dir,"/results_edgeR.tsv",sep=""), sep="\t", quote=FALSE,row.names = F, col.names=T)


