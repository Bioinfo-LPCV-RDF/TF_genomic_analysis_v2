#!/usr/bin/env Rscript
#
# R Script: Zacing_2TSS 
# =================
#
# DESCRIPTION
# This R script compute Z-scores enrichement statistic (shift from the median) 
# for spacings configurations found under peaks relative to TSS
# 
# RELATED TO compute_space
# 
# ARGUMENTS
# The script accepts the following arguments:
# - `-f, --infile`: Input file containing genomic data.
# - `-m, --mat_type`: Matrix type, either "ASYMMETRIC" or "SYMMETRIC" (default: "ASYMMETRIC").
# - `-th_tf, --thresholds_tf`: List of thresholds for TF scores (numeric).
# - `-th_tss, --thresholds_tss`: List of thresholds for TSS scores (numeric).
# - `-od, --outdir`: Output directory for results.
# - `-max`: Maximum distance to consider for spacing analysis (numeric).
# - `-c, --colors`: List of colors for plotting.
#
# FUNCTIONS
# - `foldchange`: Computes fold change for a given spacing and configuration.
# - `Zscore_calc`: Computes Z-scores for a given spacing and configuration.
# - `convert.z.score`: Converts Z-scores into p-values.
# - `myplot`: Generates Z-score vs. spacing plots with significant points labeled.
# - `volcanoPlot`: Creates volcano plots of fold change vs. p-value.
# - `myplotg`: Generates grouped Z-score vs. spacing plots for multiple thresholds.
# - `pairwise_prop_test`: Performs pairwise proportion tests between two groups.
# - `main`: Main function to compute Z-scores, generate plots, and export results for a given threshold.
#
# WORKFLOW
# 1. **Data Importation**: Reads the input file containing genomic data.
# 2. **Threshold Processing**: Iterates over user-defined thresholds to compute Z-scores and fold changes.
# 3. **Statistical Analysis**:
#    - Computes Z-scores and p-values for each spacing and configuration.
#    - Calculates pairwise proportion tests for significant spacings groups 
# 4. **Visualization**:
#    - Generates Z-score vs. spacing plots and volcano plots.
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
print("ZACING TSS: compute Zscore enrichment of motifs repeats at specific spacing, relatve to TSS")
print("===========================================")

# ------- LIBRARIES 
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library("argparse")
library(dplyr,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(cowplot,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
theme_set(theme_cowplot())
library(ggrepel,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(Cairo)

# ------- LOAD ARGUMENTS
print("----- Load arguments")
args <- commandArgs(TRUE)
parser <- ArgumentParser()

parser$add_argument("-f", "--infile")
parser$add_argument("-m", "--mat_type", default="ASYMMETRIC")
parser$add_argument('-thtf', "--thresholds_tf", type="double", nargs='+')
parser$add_argument('-thtss', "--thresholds_tss", type="double", nargs='+')
parser$add_argument('-len_thtss', "--length_thresholds_tss", type="double", nargs='+')
parser$add_argument('-len_thtf', "--length_thresholds_tf", type="double", nargs='+')
parser$add_argument("-max", type="double")
parser$add_argument("-n", "--name", type="character")

args <- parser$parse_args()

infile <- args$infile
mat_type <- args$mat_type

l_thtf <- args$length_thresholds_tf
l_thtss <- args$length_thresholds_tss
outdir <- args$name
maxs <- args$max

# extracting thresholds in lists 
thtf_list=c()
thtss_list=c()
for (i in 1:l_thtf){
	thtf_list <- append(thtf_list, as.numeric(args$thresholds_tf[i]))
}
for (i in 1:l_thtss){
	thtss_list <- append(thtss_list, as.numeric(args$thresholds_tss[i]))
}

# ------- CONFIG 
# setwd(outdir)
theme_set(theme_cowplot())
FCthreshold=0.6 # threshold for fold change to be considered significant
meth="median" #the code accepts mean and uniform but mean is bad. uniform gives similar resultat as median.

if (file.exists("Zscore_summary.txt")){
	system("rm Zscore_summary.txt")
}

# Assigning a temporary mat_type information for incoming computation (faking asymmetry even if symmetric)
mat_type_tmp<-"ASYMMETRIC"


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

myplot<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,ordmax,ordmin, colourP="grey40", colourL="black") {
# Plot Z-score vs Spacing and label significant points

	D<-data.frame(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,colourP, colourL)
		
	D$tmp<-NA
	D$tmp[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]
	#D$tmp[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05] #here we keep the negative FC..
	# with the median method I could aslo filter by Zscore> 3.5 as recommended by Iglewicz, B. and Hoaglin, D.C. (2010)
	
	res<-D$tmp[!is.na(D$tmp)]
	write(res,file="Zscore_summary.txt",append=TRUE,ncolumns=max(SPA)*2)
	
	p<-ggplot(D, aes(x=SPA, y=ZCO, label=tmp)) +
       geom_point(colour=D$colourP, size=2.5) + geom_line(colour=D$colourL) +
	   geom_text_repel() +
	   ylim(ordmin,ordmax) + 
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}

# ------------------------------------------------------------------------------

volcanoPlot<-function(SPA,FC,PVAL, ADJUSTED_PVAL) {
	# Volcano plot of Fold Change vs p-value and label significant points
	
	D<-data.frame(SPA, FC, PVAL, ADJUSTED_PVAL)
	D$FC <- ifelse(is.finite(D$FC), 0, D$FC); D$FC[is.na(D$FC)]<-0
	D$ADJUSTED_PVAL <- ifelse(is.finite(D$ADJUSTED_PVAL), 1, D$ADJUSTED_PVAL); D$ADJUSTED_PVAL[is.na(D$ADJUSTED_PVAL)]<-1
	
	D$tmp<-NA
	D$tmp[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]
	
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

pairwise_prop_test <- function(count_df, column1, group1, column2, group2, name, test_by_line=FALSE) {
  
  if (!test_by_line) { #if count are not already computed (classic df)
    tot_count <- nrow(count_df) # total population
    group1_df <- count_df[count_df[[column1]] == group1, ] # subset for group 1
    group1_count <- nrow(group1_df) # count of group 1
    group2_df <- count_df[count_df[[column2]] == group2, ] # subset for group 2
    group2_count <- nrow(group2_df) # count of group 2
    
    # proportion test
    prop_test <- prop.test(x = c(group1_count, group2_count), n = c(tot_count, tot_count))
    
    # store results
    assign(paste0("prop_test_", name), prop_test)
    assign(paste0("prop_group1_", name), as.numeric(prop_test$estimate[1]))
    assign(paste0("prop_group2_", name), as.numeric(prop_test$estimate[2]))
    assign(paste0("p.value_", name), prop_test$p.value)
    
    # create a data frame for export
	count_df_export <- data.frame(
		variable = c(paste0("prop_",group1), paste0("prop_",group2), "p.value"),
		value=c(get(paste0("prop_group1_", name)), get(paste0("prop_group2_", name)), get(paste0("p.value_", name)))
		)

    write.table(count_df_export, file=paste0("results_prop_", name,".tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  } else { #if counts are alreay computed for every lines
    
    count_df[is.na(count_df)] <- 0 # replace NA by 0 for test computation
    count_df <- count_df[rowSums(count_df[,c(3:4)]) > 0,] # remove lines without any conformations

    # proportion test for each line 
	p_value<-c()
	prop_1<-c()
	prop_2<-c()

	for (i in 1:nrow(count_df)){
	    group1_count <- count_df[[column1]][i]
        group2_count <- count_df[[column2]][i]
        tot_count <- (group1_count + group2_count)

		prop_test <- prop.test(x = c(group1_count, group2_count), n = c(tot_count, tot_count))

		assign(paste0("prop_test_", name), prop_test)
		assign(paste0("prop_", column1), as.numeric(prop_test$estimate[1]))
		assign(paste0("prop_", column2), as.numeric(prop_test$estimate[2]))
		assign(paste0("p.value_", name), prop_test$p.value)

		p_value <- append(p_value, get(paste0("p.value_", name)))
		p_value_adj <- p.adjust(p_value, method = "BH")
		prop_1 <- append(prop_1, get(paste0("prop_", column1)))
		prop_2 <- append(prop_2, get(paste0("prop_", column2)))
	}
	count_df <- cbind(count_df, p_value_adj, prop_1, prop_2)
	write.table(count_df, file=paste0("counts_",name,"_tss.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
  }

# ------------------------------------------------------------------------------

myplotg<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,ordmax,ordmin, Thresholds="black") {
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
# 	print((D))
	p<-ggplot(D,aes(x=SPA,y=ZCO, color=Thresholds,label=tmp)) +
       geom_point(size=2.5) + geom_line() +
	   geom_text_repel() +
	   ylim(ordmin,ordmax) + 
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}

# ==============================================================================

# ----- MAIN FUNCTION
main<-function(DATA, SCO, SCO2, F) {
# Main function to compute Zscore, plot and write results for a given threshold (SCO = SCORE) and a threshold name (F)
# EX: The first threshold (F1) is 1.5, the second (F2) is 2.0 and the third (F3) is 
# DATA is the input dataframe computed by bin/get_interdistances.py script in the compute_space function

	print(paste("----- Processing threshold:", SCO, "(", F, ")"))

	write("- - - - - - - - - - - - - - - - -",file="Zscore_summary.txt",append=TRUE)
	txt=paste(F, SCO, SCO2, sep="_")
	write(txt,file="Zscore_summary.txt",append=TRUE)

	filtDF<-DATA%>%filter(Score1>=SCO & Score2>=SCO2) 
# 	print(head(filtDF))

	minSpa=min(as.numeric(filtDF$Space))
	maxSpa=max(as.numeric(filtDF$Space))
	# print(minSpa)
	# print(maxSpa)
	Spacings<-seq(minSpa,maxSpa)
	outfig=paste("Zscore_", F, ".pdf", sep="")
	outtable=paste("Zscore_stats_", F, ".tsv", sep="")

	if (mat_type_tmp=="ASYMMETRIC") {
		
		FC_ER<-sapply(Spacings, foldchange, DF=filtDF, CONF="ER", S=Spacings)
		FC_IR<-sapply(Spacings, foldchange, DF=filtDF, CONF="IR", S=Spacings)
		FC_DR12<-sapply(Spacings, foldchange, DF=filtDF, CONF="DR12", S=Spacings)
		FC_DR21<-sapply(Spacings, foldchange, DF=filtDF, CONF="DR21", S=Spacings)
# 		print(head(FC_ER))
		Z_ER<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="ER", S=Spacings)
		Z_IR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="IR", S=Spacings)
		Z_DR12<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="DR12", S=Spacings)
		Z_DR21<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="DR21", S=Spacings)
# 		print(head(Z_ER))
		P_ER<-sapply(Z_ER, convert.z.score)
		P_IR<-sapply(Z_IR, convert.z.score)
		P_DR12<-sapply(Z_DR12, convert.z.score)
		P_DR21<-sapply(Z_DR21, convert.z.score)
		
		mat<-cbind(P_ER, P_IR, P_DR12, P_DR21)
		mat2<-matrix(p.adjust(as.vector(mat), method="bonferroni"),ncol=4)
		adjustedP_ER<-mat2[,1]
		adjustedP_IR<-mat2[,2]
		adjustedP_DR12<-mat2[,3]
		adjustedP_DR21<-mat2[,4]
		
		ordmax=max(max(Z_ER),max(Z_IR),max(Z_DR12),max(Z_DR21),3,na.rm=TRUE)
		ordmin=min(min(Z_ER),min(Z_IR),min(Z_DR12),min(Z_DR21),0,na.rm=TRUE)
		# print(c(ordmin,ordmax))
		
		#ER
		if (sum(is.na(mat2[,1]))) {
		p1ER<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2ER<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("ER",file="Zscore_summary.txt",append=TRUE)
		p1ER<-myplot(Spacings, Z_ER, P_ER, adjustedP_ER, FC_ER,ordmax,ordmin)
		p2ER<-volcanoPlot(Spacings, FC_ER, P_ER, adjustedP_ER)
		}
		#IR
		if (sum(is.na(mat2[,2]))) {
		p1IR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2IR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("IR",file="Zscore_summary.txt",append=TRUE)
		p1IR<-myplot(Spacings, Z_IR, P_IR, adjustedP_IR, FC_IR,ordmax,ordmin)
		p2IR<-volcanoPlot(Spacings, FC_IR, P_IR, adjustedP_IR)
		}
		#DR
		if (sum(is.na(mat2[,3]))) {
		p1DR12<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2DR12<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("DRab",file="Zscore_summary.txt",append=TRUE)
		p1DR12<-myplot(Spacings, Z_DR12, P_DR12, adjustedP_DR12, FC_DR12,ordmax,ordmin)
		p2DR12<-volcanoPlot(Spacings, FC_DR12, P_DR12, adjustedP_DR12)
		}
		#DR
		if (sum(is.na(mat2[,4]))) {
		p1DR21<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2DR21<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("DRba",file="Zscore_summary.txt",append=TRUE)
		p1DR21<-myplot(Spacings, Z_DR21, P_DR21, adjustedP_DR21, FC_DR21,ordmax,ordmin)
		p2DR21<-volcanoPlot(Spacings, FC_DR21, P_DR21, adjustedP_DR21)
		}
		
		pdf(outfig)
		print(plot_grid(p1ER, p2ER, p1IR, p2IR, p1DR12, p2DR12, p1DR21, p2DR21, ncol=2, nrow=4, labels = c('ER', '', 'IR', '', 'DRab', '', 'DRba', ''), label_size = 12))
# 		print(plot_grid(p1ER, p2ER, p1IR, p2IR, p1DR, p2DR, ncol=2, nrow=3, labels = c('ER', SCO, 'IR', '', 'DR', ''), label_size = 12))
# 		print(ggdraw()+draw_plot(p1ER, 0.0, 0.667,0.5, 0.333)+draw_plot(p1IR, 0.0, 0.334,0.5, 0.333)+draw_plot(p1DR, 0.0, 0.0,0.5, 0.333)+draw_plot(p2ER, 0.5, 0.667,0.5, 0.333)+draw_plot(p2IR, 0.5, 0.334,0.5, 0.333)+draw_plot(p2DR, 0.5, 0.0,0.5, 0.333)+draw_plot_label(c('ER', 'IR', 'DR', paste("th:",SCO,sep=" ")),c(0.0, 0.0, 0.0, 0.44),c(0.99, 0.67, 0.34, 0.99),size=12))
		dev.off()
# 		+draw_plot_label(c('ER','IR','DR',paste('th:',SCO,sep=" "),c(0.0,0.0,0.0,0.48),c(0.98,0.65,0.0,0.98),size=12))
		fi<-cbind(Spacings, FC_ER,Z_ER,P_ER,adjustedP_ER, FC_IR,Z_IR,P_IR,adjustedP_IR, FC_DR12,Z_DR12,P_DR12,adjustedP_DR12,FC_DR21,Z_DR21,P_DR21,adjustedP_DR21)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)
		

		outfig2=paste("Zscore_", F, ".png", sep="")
		cmd=paste("convert -density 144", outfig, outfig2, sep=" ")
		system(cmd)
		return(list(data.frame("Spacings"=Spacings, "FC_ER"=FC_ER,"Z_ER"=Z_ER,"P_ER"=P_ER,"Padj_ER"=adjustedP_ER, "FC_IR"=FC_IR,"Z_IR"=Z_IR,"P_IR"=P_IR,"Padj_IR"=adjustedP_IR, "FC_DR12"=FC_DR12,"Z_DR12"=Z_DR12,"P_DR12"=P_DR12,"Padj_DR12"=adjustedP_DR12, "FC_DR21"=FC_DR21,"Z_DR21"=Z_DR21,"P_DR21"=P_DR21,"Padj_DR21"=adjustedP_DR21),ordmax,ordmin))
	
	} else if (mat_type_tmp=="SYMMETRIC") {
		FC_AR<-sapply(Spacings, foldchange, DF=filtDF, CONF="AR",Spacings)
		Z_AR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="AR",Spacings)
		P_AR<-sapply(Z_AR, convert.z.score)
		adjustedP_AR<-p.adjust(P_AR, method = "bonferroni", n = length(P_AR))
	
		ordmax=max(max(Z_AR),3,na.rm=TRUE)
		ordmin=min(min(Z_AR),0,na.rm=TRUE)

		p1<-myplot(Spacings, Z_AR, P_AR, adjustedP_AR, FC_AR,ordmax,ordmin)
		p2<-volcanoPlot(Spacings, FC_AR, P_AR, adjustedP_AR)
	
		pdf(outfig, height=4, width=12)
		print(plot_grid(p1, p2, ncol=2, labels = c('A', 'B'), label_size = 12))
		dev.off()
		
		fi<-cbind(Spacings, FC_AR,Z_AR,P_AR,adjustedP_AR)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)

		outfig2=paste("Zscore_", F, ".png", sep="")
		cmd=paste("convert -density 144", outfig, outfig2, sep=" ")
		system(cmd)
		return(list(data.frame("Spacings"=Spacings, "FC_AR"=FC_AR,"Z_AR"=Z_AR,"P_AR"=P_AR,"Padj_AR"=adjustedP_AR),ordmax,ordmin))
	} # AR: All Repeats
}

#infile look like this:
#Peak	Spacing	Score1	Score2	matricePosition1	matricePosition2	correctedPosition1	correctedPosition2
#chr1:19257041-19257441	DR_6	-4.49980967033	-5.19295685089	6	23	17	23
#chr1:19257041-19257441	DR_11	-4.49980967033	-5.29831736655	6	28	17	28

d<-read.table(file=infile, header=TRUE)

setwd(outdir)

# print(nrow(d))

# print(head(d))

# d2 <- d %>% 
# 	separate(Peak, c("tmp", "End"), sep="-") %>%
# 	separate(tmp, c("chr", "Start"), sep=":") %>%
# 	separate(Spacing,c("Conf", "Space"), sep="_")%>%
# 	mutate(Size=as.integer(End)-as.integer(Start)+1)

tx="Listing of over or under represented Spacings for selected scores thresholds"
write(tx,file="Zscore_summary.txt",append=TRUE)

# Definitions before for loop
DF <- NULL
list_F <- c()
returnF <- c()
ordmax=c()
ordmin=c()

# If only 1 threshold was chosen for TSS from spacing_TF2TSS
if (l_thtss < l_thtf) {
	for (i in 1:l_thtf){
		thtss_list <- append(thtss_list, thtss_list[1])
	}
}

# print(thtss_list)

for (i in 1:length(thtf_list)) {
  # create a returnF variable for every threshold and pass them into the main function (to give as 3 dfs)
  returnF_var <- paste("returnF", i, sep="")
  assign(returnF_var, main(d, thtf_list[[i]], thtss_list[[i]], paste("F", i, sep="")))

  # create a list_F variable for every threshold to extract the first df in correspondng returnF variable
  list_F_var <- paste("list_F", i, sep="")
  assign(list_F_var, get(returnF_var)[[1]])

  # adding a Thresholds column in df
  list_F_var_df <- get(list_F_var)  # Extract the data frame
  list_F_var_df$Thresholds <- paste("TF(a): ",round(thtf_list[[i]], 2), "; TSS(b): ",round(thtss_list[[i]], 2))
  assign(list_F_var, list_F_var_df)

  # chosing ordmax and min (2nd and 3rd dfs in returnF variables)
  ordmax_F_var <- paste("ordmax_F", i, sep="")
  ordmax_F_var <- get(returnF_var)[[2]]
  ordmax <- append(ordmax, ordmax_F_var)
  ordmin_F_var <- paste("ordmin_F", i, sep="")
  ordmin_F_var <- get(returnF_var)[[3]]
  ordmin <- append(ordmin, ordmin_F_var)
	}

# chosing maximum and minimum ord of the different threshold computated dfs
ordmax <- max(ordmax)
ordmin <- min(ordmin)

# merging every thing into a new DF
DF <- do.call(rbind, mget(paste0("list_F", 1:length(thtf_list))))
DF$Thresholds <- as.factor(DF$Thresholds)

# print(DF)
if (mat_type_tmp=="ASYMMETRIC") {
#ER
p1ER<-myplotg(DF$Spacings, DF$Z_ER, DF$P_ER, DF$Padj_ER, DF$FC_ER,ordmax,ordmin, DF$Thresholds)
print("ER plot done")
#IR
p1IR<-myplotg(DF$Spacings, DF$Z_IR, DF$P_IR, DF$Padj_IR, DF$FC_IR,ordmax,ordmin, DF$Thresholds)
print("IR plot done")
#DR
p1DR12<-myplotg(DF$Spacings, DF$Z_DR12, DF$P_DR12, DF$Padj_DR12, DF$FC_DR12,ordmax,ordmin, DF$Thresholds)
print("DRab plot done")
#DR
p1DR21<-myplotg(DF$Spacings, DF$Z_DR21, DF$P_DR21, DF$Padj_DR21, DF$FC_DR21,ordmax,ordmin, DF$Thresholds)
print("DRba plot done")

pdf("Zscore_allF.pdf")
print(plot_grid(p1ER, p1IR, p1DR12, p1DR21, ncol=1, nrow=4, labels = c('ER', 'IR', 'DRab', 'DRba', ''), label_size = 12))
dev.off()

} else if (mat_type_tmp=="SYMMETRIC") {

p1AR<-myplotg(DF$Spacings, DF$Z_AR, DF$P_AR, DF$Padj_AR, DF$FC_AR,ordmax,ordmin, DF$Thresholds)
print("AR plot done")

pdf("Zscore_allF.pdf")
print(plot_grid(p1AR, ncol=1, nrow=4, labels = c('AR'), label_size = 12))
dev.off()

}

cmd="convert -density 144 Zscore_allF.pdf Zscore_allF.png"
system(cmd)

# df1 <- data.frame(matrix(unlist(F1_list), nrow=length(F1_list), byrow=TRUE))
# df2 <- data.frame(matrix(unlist(F2_list), nrow=length(F2_list), byrow=TRUE))
# df3 <- data.frame(matrix(unlist(F3_list), nrow=length(F3_list), byrow=TRUE))
# create a summary_plot with all scores

# extract candidates that are detected in two of three Fscore selections (note uniq -d means at least two)
if (mat_type_tmp=="ASYMMETRIC") {
system("cat Zscore_summary.txt | grep ER -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/ER/g' > candidates.txt")
system("cat Zscore_summary.txt | grep IR -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/IR/g' >> candidates.txt")
system("cat Zscore_summary.txt | grep DRab -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/DRab/g' >> candidates.txt")
system("cat Zscore_summary.txt | grep DRba -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/DRba/g' >> candidates.txt")
} else if (mat_type_tmp=="SYMMETRIC") {
system("cat Zscore_summary.txt | grep '[1-9]' | tr ' ' '\n' | sort |
uniq -d > candidates.txt")
}
# warnings()


# ----------------------- Visualisation TSS spacings/orientation -------------------------

# This is the Spacing_TFBSTSS visual representation and statistics script. Informations are 
# extracted from candidates.txt, created previously. 

# 1 - Extracting scpacing values from candidates.txt

if (file.size("candidates.txt")!=0){ #if candidates file is not empty 

	list_candidates <- read.table(file="candidates.txt")
	list_candidates <- list_candidates$V1
	print("Candidates liste:")
	print(list_candidates)

	# 2 - Getting data for the TSS graph 

	# Initialisation (creating lists which will be df columns, storing each candidates informations)

	candidate_list <- c() #listing the differents significant candidates -> from candidates.txt
	spacing_list <- c() #spacings conformations (DRab, DRba, IR, ER) -> from candidates.txt
	distance_list <- c() #distances of TFBS from TSS -> from candidates.txt
	Zvalue_list <- c() #values of Zscores computation for every thresholds -> from Zscores_stats_F
	Psig_list <- c() #values of Zscores computation for every thresholds -> from Zscores_stats_F
	orientation_list <- c() #either TRUE if TFBS is in the same direction of TSS or FALSE 
	position_list <- c() #either "before" if TFBS is before TSS or "after"

	# Extracting data from candidates list:
	for (i in 1:length(list_candidates)) {

		candidate_n <- list_candidates[i] #extracting each candidate

		spacing_var <- paste("spacing_n", i, sep = "")
		assign(spacing_var, gsub("[[:digit:]]", "", candidate_n)) #extracting from each candidate its corresponding spacing

		distance_var <- paste("distance_n", i, sep = "")
		assign(distance_var, gsub("[^0-9.-]", "", candidate_n)) #extracting from each candidate its corresponding distance

		
		# Searching Zscore table from outdir
		for (j in 1:length(thtf_list)) {
			Zscore_stats_F_var <- paste("Zscore_stats_F", j, ".tsv", sep = "")
			Score_F_var <- paste("Score_F", j, sep = "")
			assign(Score_F_var, read.table(file = Zscore_stats_F_var, sep = "\t", header = TRUE))
			
			# Extracting row for each candidate corresponding distance
			subset_score <- get(Score_F_var)[get(Score_F_var)$Spacings == as.numeric(get(distance_var)), ]

			# Extracting candidates informations depending on values from Zscore_stats_F columns files
			Zvalue_var <- subset_score[[if (get(spacing_var) == "DRab") 11 else if (get(spacing_var) == "DRba") 15 else if (get(spacing_var) == "IR") 7 else if (get(spacing_var) == "AR") 3 else 3]]
			Psig_var <- subset_score[[if (get(spacing_var) == "DRab") 13 else if (get(spacing_var) == "DRba") 17 else if (get(spacing_var) == "IR") 9 else if (get(spacing_var) == "AR") 5 else 5]]
			
			orientation_var <- ifelse(get(spacing_var) %in% c("DRab","DRba"), TRUE, FALSE)
			position_var <- ifelse(get(spacing_var) %in% c("DRab", "ER"), "before", "after")

			# Assigning variables to lists
			candidate_list <- append(candidate_list, rep(candidate_n, nrow(subset_score)))
			spacing_list <- append(spacing_list, rep(get(spacing_var), nrow(subset_score)))
			distance_list <- append(distance_list, rep(get(distance_var), nrow(subset_score)))
			Zvalue_list <- append(Zvalue_list, Zvalue_var)
			Psig_list <- append(Psig_list, Psig_var)
			orientation_list <- append(orientation_list, rep(orientation_var, nrow(subset_score)))
			position_list <- append(position_list, rep(position_var, nrow(subset_score)))
		}
	}

	# 3 - Merging all extracted data into a df 

	df_TSS <- data.frame(
	candidates = candidate_list,
	distance = distance_list,
	spacing = spacing_list,
	orientation = orientation_list,
	position = position_list,
	Zvalue = Zvalue_list,
	Psig = Psig_list
	)

	# computing mean Zvalue 
	df_TSS$Zvalue[which(df_TSS$Zvalue == "-Inf")] <- NA #explicit NA
	df_TSS$Zvalue[which(df_TSS$Zvalue == "Inf")] <- NA #explicit NA
	mean_zvalue <- aggregate(Zvalue ~ candidates, data = df_TSS, FUN = mean) #compute mean
	df_TSS <- merge(df_TSS, mean_zvalue, by ="candidates") #mean zvalues is Zvalue.y column

	#  if mean Zvalue is negative (non enriched), remove it from the df
	df_TSS <- df_TSS %>%
	filter(Zvalue.y > 0)

	# If TFBS is enriched before TSS (position is "before"),  it should be represented with negative distance 
		df_TSS <- df_TSS %>%
		mutate(distance = ifelse(as.character(position)=="before", as.numeric(as.character(df_TSS$distance)) * (-1), as.numeric(as.character(df_TSS$distance))))
	
	# But if matrix is SYMMETRIC, configuration should be converted into AR (All repeats)
	if (mat_type=="SYMMETRIC"){
		df_TSS <- df_TSS %>%
			distinct(distance, .keep_all = TRUE) %>% #distance needed only once
			mutate(spacing="AR", orientation="NONE")
	}

	# reduce the graph to only get one occurence for each candidate (only used for labels in geom_text_repel)
	df_TSS_unique <- df_TSS[,-c(6:7)] #remove Zvalue.x and Psig
	df_TSS_unique <- df_TSS_unique[!duplicated(df_TSS_unique), ] #remove duplicated rows

	# exporting df into a table (used for larger analysis, ex : TSS analysis per family)
	print("Exporting candidates informations to candidates_for_graph.tsv")
	# Example : 
	# candidates	distance	spacing	orientation	position	Zvalue.y
	# 1	-41	DRab	TRUE	before	4.72143531633617
	print(df_TSS_unique)
	write.table(df_TSS_unique, file="candidates_for_graph.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


	# 4 - Statistics and data frame manipulation before visualisation

	# global enrichment : test 2 groups proportion at a global level (usually within the global 1000 bp)
	if (mat_type=="ASYMMETRIC"){
		pairwise_prop_test(df_TSS_unique, "orientation", "FALSE", "orientation", "TRUE", "orientation", test_by_line=FALSE) #function defined earlier
		print("Global TFBS enrichment results save as results_prop_orientation.tsv")

		# proportion test within distance ranges (sliding windows of 10 bp)
		count_before_df <- read.table(file="counts_before_tss.tsv", header=TRUE, sep="\t")
		count_after_df <- read.table(file="counts_after_tss.tsv", header=TRUE, sep="\t")
		pairwise_prop_test(count_before_df, column1="count_DRab", column2="count_ER", name="before", test_by_line=TRUE)
		pairwise_prop_test(count_after_df, column1="count_DRba", column2="count_IR", name="after", test_by_line=TRUE)
		print("Orientation enrichment results save as counts_after_tss.tsv and counts_before_tss.tsv")

	}

	# if SYMMETRIC matrix = only global position (before/after TSS) propportionaly test is relevant 
	pairwise_prop_test(df_TSS_unique, "position", "before", "position", "after", "position", test_by_line=FALSE)
	print("Global TFBS enrichment results save as results_prop_position.tsv")


	#  3 - Plot

	p <- ggplot(df_TSS, aes(x = factor(orientation), y = distance, shape = orientation, label = distance, color=Zvalue.y, alpha=Zvalue.y)) +

		geom_segment(aes(x = 0, y = 0, xend = length(unique(df_TSS_unique$orientation)) +0.5, yend = 0), col = "black") +
		geom_segment(aes(x = length(unique(df_TSS_unique$orientation)) +0.5, y = 0, xend = length(unique(df_TSS_unique$orientation)) + 0.5, yend = 90),
					col = "black",
					arrow = arrow(length = unit(0.5, "cm"))) +
		
		geom_point(aes(x = factor(orientation), y = distance), size = 15) +
		scale_color_gradient(low = alpha("red", 0), high = "red", limits = c(0,NA)) +
		coord_flip() +
		scale_shape_manual(values = c("TRUE" = "►", "FALSE" = "◄", "NONE" = "◆"), labels = c("Same", "Opposite", "Palindromic")) +
		
		geom_text(label = "TSS", y = 150, x = length(unique(df_TSS_unique$orientation)) +0.5, col = "black") +
		geom_text_repel(data = df_TSS_unique, size = 3) +
		ggtitle("TFBS preferential distance and orientation from TSS\n") +
		ylim(-maxs, maxs) +
		labs(col = "Z-scores mean", shape = "Orientation from TSS") +
		scale_x_discrete(labels = NULL, breaks = NULL) +
		xlab("") +
		ylab("Distance from TSS (bp)") +
		theme(axis.text.x = element_text(size = 15))+
		guides(alpha="none")

	
	ggsave("TSS_spacing.pdf", p, width = 15, height = 4, device = cairo_pdf)
	ggsave("TSS_spacing.png", p, width = 15, height = 4, type = "cairo")
	print("TSS spacings plotted at TSS_spacing.pdf and TSS_spacing.png")

} else { #if candidates file is empty
	print("There is no candidates for this spacing : try to increase $filter argument if not 1000 bp")
	quit()
}




