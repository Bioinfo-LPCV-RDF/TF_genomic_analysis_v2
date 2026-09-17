#!/usr/bin/env Rscript
#
# R Script: Zacing 
# =================
#
# DESCRIPTION
# This R script compute Z-scores enrichement statistic (shift from the median) 
# for spacings configurations found under peaks
# 
# RELATED TO compute_space
# 
# ARGUMENTS
# The script accepts the following arguments:
# - `-f, --infile`: Input file containing genomic data.
# - `-m, --mat_type`: Matrix type, either "ASYMMETRIC" or "SYMMETRIC" (default: "ASYMMETRIC").
# - `-th, --thresholds`: List of thresholds for filtering scores (type: double).
# - `-od, --outdir`: Output directory for results.
# - `--maxy`: Maximum y-axis value for plots (optional).
# - `--miny`: Minimum y-axis value for plots (optional).
# - `-c, --colors`: List of colors for plotting.
# - `-random`: Number of random iterations for empirical p-value computation (default: 0).
# - `-correction`: Method for p-value adjustment (e.g., "bonferroni").
#
# FUNCTIONS
# - `foldchange`: Computes fold change for a given spacing and configuration.
# - `Zscore_calc`: Computes Z-scores for a given spacing and configuration.
# - `convert.z.score`: Converts Z-scores into p-values.
# - `myplot`: Generates Z-score vs. spacing plots with significant points labeled.
# - `volcanoPlot`: Creates volcano plots of fold change vs. p-value.
# - `myplotg`: Generates grouped Z-score vs. spacing plots for multiple thresholds.
# - `main`: Main function to compute Z-scores, generate plots, and export results for a given threshold.
#
# WORKFLOW
# 1. **Data Importation**: Reads the input file containing genomic data.
# 2. **Threshold Processing**: Iterates over user-defined thresholds to compute Z-scores and fold changes.
# 3. **Statistical Analysis**:
#    - Computes theoretical or empirical p-values based on the `-random` argument.
#    - Adjusts p-values using the specified correction method.
# 4. **Visualization**:
#    - Generates Z-score vs. spacing plots and volcano plots.
#    - Adds confidence intervals if empirical p-values are computed.
# 5. **Results Exportation**:
#    - Saves results as tables and plots in the specified output directory.
#    - Identifies significant spacings and writes them to summary files.
#
# OUTPUT
# - `Zscore_allF.pdf`: Combined plots for all thresholds.
# - `Zscore_summary.txt`: Summary of significant spacings.
# - `candidates.txt`: List of candidate spacings detected in multiple thresholds.
#
# NOTES
# - The script supports both symmetric and asymmetric matrix types.
#

# ----- LEXICAL NOTES
# BP: Distance between repeats (in bp)
# CONF: Configuration of the repeats (ER, IR, DR, AR)
# SPA: SPACING (a CONFIGURATION and a distance in BP association)
# ER: Everted Repeats (<--  -->)
# DR: Direct Repeats (-->  --> or <--  <--)
# IR: Inverted Repeats (-->  <--)
# AR: All Repeats (<-->  <-->, for symetric motif with no inherent orientation)
# FC: Fold Change
# Z: Z-score computed from FC

print("===========================================")
print("ZACING: compute Zscore enrichment of motifs repeats at specific spacing")
print("===========================================")

# ------- LIBRARIES 
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library("argparse")
library(dplyr,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(cowplot,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggrepel,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)

# ------- LOAD ARGUMENTS
print("----- Load arguments")
args <- commandArgs(TRUE)
parser <- ArgumentParser()

parser$add_argument("-f", "--infile")
parser$add_argument("-m", "--mat_type", default="ASYMMETRIC")
parser$add_argument('-th', "--thresholds", type="double", nargs='+')
parser$add_argument("-od", "--outdir")
parser$add_argument("--maxy", type="double")
parser$add_argument("--miny", type="double")
parser$add_argument('-c', "--colors", nargs='+')

args <- parser$parse_args()
infile <- args$infile
mat_type <- args$mat_type
ths_list <- args$thresholds
outdir <- args$outdir
maxy <- args$maxy
miny <- args$miny
colors <- args$colors

# print(args)

# ------- CONFIG 
theme_set(theme_cowplot())
FCthreshold=0.6 # threshold for fold change to be considered significant
meth="median" #the code accepts mean and uniform but mean is bad. uniform gives similar resultat as median.

if (file.exists("Zscore_summary.txt")){
	system("rm Zscore_summary.txt")
}


# ----- FUNCTIONS
# ------------------------------------------------------------------------------
foldchange<-function(BP, CONF, DF, S) {
# Compute fold change for a given spacing (BP) and configuration (CONF)

	if (CONF=="AR"){
		DF<-DF
	} else {
		DF<-DF %>% filter(Conf==CONF)
	}
	x<-nrow(DF%>%filter(Space==BP))
	Scount<-table(DF$Space) # this is to compute mean and median
	
	if (meth=="uniform") {
		population=nrow(DF)
		u=population*(1/length(S)) # same as population * (max - min)/2
	}
	else if (meth=="mean") {u=mean(Scount)}
	else if (meth=="median") {u=median(Scount)}
	
	FC=x/u

	return(FC)	
}

# ------------------------------------------------------------------------------

Zscore_calc<-function(BP, CONF, DF, S, meth="median") {
# Compute Z-score for a given spacing (BP) and configuration (CONF)
#formula adpated from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142402/

	if (CONF=="AR"){
		DF<-DF
	} else {
		DF<-DF %>% filter(Conf==CONF)
	}
	
	x<-nrow(DF%>%filter(Space==BP))
	Scount<-table(DF$Space)
	
	if (meth=="uniform") {
		population=nrow(DF)
		u=population*(1/length(S)) # same as population * (max - min)/2

		P=1/length(S)
		sigma=sqrt(population*P*(1-P))

		#Calculations variations (not exactly same results)
		#variance=((max(S)-min(S))^2)/12
		#sigma=sqrt(variance)

	}
	
	else if (meth=="mean") {
		# compute the real u from mean or median 
		u=mean(Scount)
		sigma=sd(Scount)
	}
	else if (meth=="median") {
		# mad computes the median absolute deviation
		u=median(Scount)
		sigma=mad(Scount) 
	}
	
	Z=(x-u)/sigma	
	return(Z)
}

# ------------------------------------------------------------------------------

convert.z.score<-function(z, one.sided=NULL) {
# Convert a Z-score into a p-value
#https://www.biostars.org/p/17227/

    if(is.null(one.sided)) {
        pval = pnorm(-abs(z));
        pval = 2 * pval
    } else if(one.sided=="-") {
        pval = pnorm(z);
    } else {
        pval = pnorm(-z);                                                                                 
    }
    return(pval);
}  

# ------------------------------------------------------------------------------

myplot<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,ordmax,ordmin) {
# Plot Z-score vs Spacing and label significant points

	D<-data.frame(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC)
		
	D$tmp<-NA
	D$tmp[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]
	#D$tmp[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05] #here we keep the negative FC..
	# with the median method I could aslo filter by Zscore> 3.5 as recommended by Iglewicz, B. and Hoaglin, D.C. (2010)
	
	res<-D$tmp[!is.na(D$tmp)]
	write(res,file="Zscore_summary.txt",append=TRUE)
	
	p<-ggplot(D, aes(x=SPA, y=ZCO, label=tmp)) +
       geom_point(colour="grey40", size=2.5) + geom_line() +
	   geom_text_repel() +
	   ylim(ordmin,ordmax) + 
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}

# ------------------------------------------------------------------------------

volcanoPlot<-function(SPA,FC,PVAL, ADJUSTED_PVAL) {
# Volcano plot of Fold Change vs p-value and label significant points
	
	D<-data.frame(SPA, FC, PVAL, ADJUSTED_PVAL)
	
	D$tmp<-NA
	D$tmp[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]
	#D$tmp[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]
	myCol="darkblue"
	
	
	p<-ggplot(D,aes(x=log2(FC), y=-log10(ADJUSTED_PVAL), label=tmp)) + 
		geom_point() + 
		geom_text_repel() +
		xlab("log2(FC)") + ylab("-log10(P)") +
		geom_vline(xintercept=c(-FCthreshold, FCthreshold), col=myCol, linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), col=myCol, linetype="dashed")
	return(p)
}

# ------------------------------------------------------------------------------

myplotg<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,ordmax,ordmin, Thresholds="black", colors) {
# Plot Z-score vs Spacing and label significant points, with multiple thresholds
	D<-data.frame(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC, Thresholds)
		
	D$ADJUSTED_PVAL <- ifelse(is.finite(D$ZCO), D$ADJUSTED_PVAL,1)
	D$FC <- ifelse(is.finite(D$ZCO), D$FC,1)
	D$ZCO <- ifelse(is.finite(D$ZCO), D$ZCO,0); D$ZCO[is.na(D$ZCO)]<-0
	D$FC <- ifelse(is.finite(D$FC), D$FC,0); D$FC[is.na(D$FC)]<-0
	D$ADJUSTED_PVAL <- ifelse(is.finite(D$ADJUSTED_PVAL), D$ADJUSTED_PVAL,1); D$ADJUSTED_PVAL[is.na(D$ADJUSTED_PVAL)]<-1
	
	D$tmp<-NA
	D$tmp[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]
	# with the median method I could aslo filter by Zscore> 3.5 as recommended by Iglewicz, B. and Hoaglin, D.C. (2010)

	p<-ggplot(D,aes(x=SPA,y=ZCO, color=Thresholds,label=tmp)) +
       geom_point(size=2.5) + geom_line() +
	   geom_text_repel() +
	   ylim(ordmin,ordmax) + 
	   scale_color_manual(values=colors)+
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}


# ==============================================================================

# ----- MAIN FUNCTION
main<-function(DATA, SCO, F) {
# Main function to compute Zscore, plot and write results for a given threshold (SCO = SCORE) and a threshold name (F)
# EX: The first threshold (F1) is 1.5, the second (F2) is 2.0 and the third (F3) is 
# DATA is the input dataframe computed by bin/get_interdistances.py script in the compute_space function

	print(paste("----- Processing threshold:", SCO, "(", F, ")"))

	txt=paste(F, SCO, sep="_")
	write(txt,file="Zscore_summary.txt",append=TRUE)

	filtDF<-DATA%>%filter(Score1>SCO & Score2>SCO) # keep only motifs with scores above threshold

	minSpa=min(as.numeric(filtDF$Space))
	maxSpa=max(as.numeric(filtDF$Space))
	Spacings<-seq(minSpa,maxSpa)
	outfig=paste("Zscore_", F, ".pdf", sep="")
	outtable=paste("Zscore_stats_", F, ".tsv", sep="")
	Znull_storage <- list()
	
	if (mat_type=="ASYMMETRIC") { #ER, IR, DR
		print("-> Selected motif is asymmetric")
		print("-> Computing Z-scores...")
		FC_ER<-sapply(Spacings, foldchange, DF=filtDF, CONF="ER", S=Spacings)
		FC_IR<-sapply(Spacings, foldchange, DF=filtDF, CONF="IR", S=Spacings)
		FC_DR<-sapply(Spacings, foldchange, DF=filtDF, CONF="DR", S=Spacings)
		Z_ER<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="ER", S=Spacings)
		Z_IR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="IR", S=Spacings)
		Z_DR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="DR", S=Spacings)

		print("-> Computing p-values...")
		P_ER<-sapply(Z_ER, convert.z.score)
		P_IR<-sapply(Z_IR, convert.z.score)
		P_DR<-sapply(Z_DR, convert.z.score)
		mat<-cbind(P_ER, P_IR, P_DR)
		pvalues<-matrix(p.adjust(as.vector(mat), method="bonferroni"),ncol=3)
		adjustedP_ER<-pvalues[,1]
		adjustedP_IR<-pvalues[,2]
		adjustedP_DR<-pvalues[,3]

		# export table
		print("-> Exporting table...")
		fi<-cbind(Spacings, FC_ER,Z_ER,P_ER,adjustedP_ER, FC_IR,Z_IR,P_IR,adjustedP_IR, FC_DR,Z_DR,P_DR,adjustedP_DR)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)

		# determine max and min for y axis
		if (maxy==0){
			ordmax=max(max(Z_ER),max(Z_IR),max(Z_DR),3,na.rm=TRUE)
		}else{
			ordmax=maxy
		}
		
		if(miny==0){
			ordmin=min(min(Z_ER),min(Z_IR),min(Z_DR),0,na.rm=TRUE)
		}else{
			ordmin=miny
		}

		print("-> Plotting results...")
		if (sum(is.na(pvalues[,1]))) {
			p1ER<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
			p2ER<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		} else
		{
			write("ER",file="Zscore_summary.txt",append=TRUE)
			p1ER<-myplot(Spacings, Z_ER, P_ER, adjustedP_ER, FC_ER,ordmax,ordmin)
			p2ER<-volcanoPlot(Spacings, FC_ER, P_ER, adjustedP_ER)
		}
		if (sum(is.na(pvalues[,2]))) {
			p1IR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
			p2IR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
			write("IR",file="Zscore_summary.txt",append=TRUE)
			p1IR<-myplot(Spacings, Z_IR, P_IR, adjustedP_IR, FC_IR,ordmax,ordmin)
			p2IR<-volcanoPlot(Spacings, FC_IR, P_IR, adjustedP_IR)
		}
		if (sum(is.na(pvalues[,3]))) {
			p1DR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
			p2DR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
			write("DR",file="Zscore_summary.txt",append=TRUE)
			p1DR<-myplot(Spacings, Z_DR, adjustedP_DR, adjustedP_DR, FC_DR,ordmax,ordmin)
			p2DR<-volcanoPlot(Spacings, FC_DR, adjustedP_DR, adjustedP_DR)
		}

		# export more detailled plot with volcano plots
		pdf(outfig)
		print(ggdraw()+draw_plot(p1ER, 0.0, 0.667,0.5, 0.333)+draw_plot(p1IR, 0.0, 0.334,0.5, 0.333)+draw_plot(p1DR, 0.0, 0.0,0.5, 0.333)+draw_plot(p2ER, 0.5, 0.667,0.5, 0.333)+draw_plot(p2IR, 0.5, 0.334,0.5, 0.333)+draw_plot(p2DR, 0.5, 0.0,0.5, 0.333)+draw_plot_label(c('ER', 'IR', 'DR', paste("th:",SCO,sep=" ")),c(0.0, 0.0, 0.0, 0.39),c(0.99, 0.67, 0.34, 0.99),size=12))
		dev.off()
		OUT=list(data.frame("Spacings"=Spacings, "FC_ER"=FC_ER,"Z_ER"=Z_ER,"P_ER"=P_ER,"Padj_ER"=adjustedP_ER, "FC_IR"=FC_IR,"Z_IR"=Z_IR,"P_IR"=P_IR,"Padj_IR"=adjustedP_IR, "FC_DR"=FC_DR,"Z_DR"=Z_DR,"P_DR"=P_DR,"Padj_DR"=adjustedP_DR),ordmax,ordmin, Znull_storage)

		
	} else if (mat_type=="SYMMETRIC") { # AR
		print("-> Symmetric mode")
		print("-> Computing Z-scores...")
		FC_AR <- sapply(Spacings, foldchange, DF=filtDF, CONF="AR", S=Spacings)
		Z_AR  <- sapply(Spacings, Zscore_calc, DF=filtDF, CONF="AR", S=Spacings)
		
		print("-> Computing p-values...")
		P_AR <- sapply(Z_AR, convert.z.score)
		adjustedP_AR <- p.adjust(P_AR, method = "bonferroni", n = length(P_AR))

		print("-> Exporting table...")
		fi <- cbind(Spacings, FC_AR, Z_AR, P_AR, adjustedP_AR)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names=FALSE)
		
		ordmax <- if(maxy==0) max(max(Z_AR, na.rm=TRUE),3) else maxy
		ordmin <- if(miny==0) min(min(Z_AR, na.rm=TRUE),0) else miny
		
		# Plots
		print("-> Plotting results...")
		if (sum(is.na(P_AR))) {
			p1 <- ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
			p2 <- ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		} else {
			p1 <- myplot(Spacings, Z_AR, P_AR, adjustedP_AR, FC_AR, ordmax, ordmin)
			p2 <- volcanoPlot(Spacings, FC_AR, P_AR, adjustedP_AR)
		}

		pdf(outfig, height=4, width=12)
		print(plot_grid(p1, p2, ncol=2, labels = c('A', 'B'), label_size = 12))
		dev.off()
		
		OUT <- list(data.frame("Spacings"=Spacings, "FC_AR"=FC_AR, "Z_AR"=Z_AR, "P_AR"=P_AR, "Padj_AR"=adjustedP_AR),
					ordmax, ordmin,Znull_storage)
	} 
	return(OUT)
}


# ==============================================================================


# ----- DATA IMPORTATION
d<-read.table(infile, header=TRUE) 
#infile look like this:
#Peak	Spacing	Score1	Score2	matricePosition1	matricePosition2	correctedPosition1	correctedPosition2
#chr1:19257041-19257441	DR_6	-4.49980967033	-5.19295685089	6	23	17	23
#chr1:19257041-19257441	DR_11	-4.49980967033	-5.29831736655	6	28	17	28

setwd(outdir)

pvalue_method <- "theoretical p-values"

tx=(paste0("Listing of over or under represented Spacings for selected scores thresholds (", pvalue_method, ")"))
write(tx,file="Zscore_summary.txt",append=TRUE)

# ----- FUNCTIONS CALLS
print("LAUCHING - Zscore computation")
# List definitions before for loop
DF <- NULL
list_F <- c()
returnF <- c()
ordmax=c()
ordmin=c()

for (i in 1:length(ths_list)) {
  # create a returnF variable for every threshold and pass them into the main function (to give as 3 dfs)
  returnF_var <- paste("returnF", i, sep="")
  assign(returnF_var, main(d, ths_list[[i]], paste("F", i, sep="")))

  # create a list_F variable for every threshold to extract the first df in corresponding returnF variable
  list_F_var <- paste("list_F", i, sep="")
  assign(list_F_var, get(returnF_var)[[1]])

  # adding a Thresholds column in df
  list_F_var_df <- get(list_F_var)  # Extract the data frame
  list_F_var_df$Thresholds <- round(ths_list[[i]], 2)
  assign(list_F_var, list_F_var_df)

  # chosing ordmax and min (2nd and 3rd dfs in returnF variables)
  ordmax_F_var <- paste("ordmax_F", i, sep="")
  ordmax_F_var <- get(returnF_var)[[2]]
  ordmax <- append(ordmax, ordmax_F_var)
  ordmin_F_var <- paste("ordmin_F", i, sep="")
  ordmin_F_var <- get(returnF_var)[[3]]
  ordmin <- append(ordmin, ordmin_F_var)
	}

ordmax <- max(ordmax)
ordmin <- min(ordmin)

# ----- FINAL RESULTS EXPORTATION: Plots and Tables 
print("----- Exporting final results...")
DF <- do.call(rbind, mget(paste0("list_F", 1:length(ths_list))))
DF$Thresholds <- as.factor(DF$Thresholds)

if (mat_type == "ASYMMETRIC") {
  # Plots de base pour chaque configuration
  p1ER <- myplotg(DF$Spacings, DF$Z_ER, DF$P_ER, DF$Padj_ER, DF$FC_ER, ordmax, ordmin, DF$Thresholds, colors)
  p1IR <- myplotg(DF$Spacings, DF$Z_IR, DF$P_IR, DF$Padj_IR, DF$FC_IR, ordmax, ordmin, DF$Thresholds, colors)
  p1DR <- myplotg(DF$Spacings, DF$Z_DR, DF$P_DR, DF$Padj_DR, DF$FC_DR, ordmax, ordmin, DF$Thresholds, colors)

  # Export PDF final
  pdf("Zscore_allF.pdf")
  print(plot_grid(p1ER, p1IR, p1DR, ncol=1, nrow=3, labels = c('ER', 'IR', 'DR'), label_size = 12))
  dev.off()

}  else if (mat_type=="SYMMETRIC") {
	p1AR<-myplotg(DF$Spacings, DF$Z_AR, DF$P_AR, DF$Padj_AR, DF$FC_AR,ordmax,ordmin, DF$Thresholds, colors)

	pdf("Zscore_allF.pdf")
	print(plot_grid(p1AR, ncol=1, nrow=1, labels = c('AR'), label_size = 12))
	dev.off()
} 

cmd="pdftoppm Zscore_allF.pdf Zscore_allF -png; mv Zscore_allF-1.png Zscore_allF.png"
system(cmd)

# extract candidates that are detected in two of three Fscore selections (note uniq -d means at least two)
if (mat_type=="ASYMMETRIC") {
system("cat Zscore_summary.txt | grep ER -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/ER/g' > candidates.txt")
system("cat Zscore_summary.txt | grep IR -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/IR/g' >> candidates.txt")
system("cat Zscore_summary.txt | grep DR -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/DR/g' >> candidates.txt")
} else if (mat_type=="SYMMETRIC") {
system("cat Zscore_summary.txt | grep '[0-9]' | tr ' ' '\n' | sort |
uniq -d > candidates.txt")
}

print("ALL DONE: Zscore computation finished, results in Zscore_allF.pdf and Zscore_summary.txt")