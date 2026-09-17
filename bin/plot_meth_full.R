#' Methylation Analysis and Plotting Script
#' ========================================
#' 
#' DESCRIPTION
#' This script performs analysis and visualization of methylation data 
#' in transcription factor binding sites (TFBS). It calculates correlations, generates 
#' plots, and saves results for further analysis.
#' 
#' RELATED TO cooking_methylation 
#' 
#' PROCESS
#' The script performs the following tasks:
#' 1. Loads required libraries and initializes the working environment.
#' 2. Reads input arguments and prepares the dataset for analysis.
#' 3. Defines a function `plot_and_test` to calculate correlations and generate scatter plots.
#' 4. Loads the dataset and computes additional metrics for analysis.
#' 5. Generates histograms and scatter plots for methylation data.
#' 6. Performs per-position regression analysis for methylation probabilities.
#' 7. Computes correlation matrices for methylation probabilities between positions.
#' 8. Saves results and generates output files for downstream analysis.
#' 
#' ARGUMENTS
#' - [1] Path to the working directory containing input files.
#' - [2] Motif length (integer) for the transcription factor binding site.
#' - [3] Motif symmetry ("yes" or "no").
#' - [4] Score cutoff (integer) to filter TFBS based on binding scores.
#' 
#' INPUT
#' - `table_peak_bs_meth_Zhu.txt`: Input file containing methylation data.
#' - `vectorsForFigs.RData`: File to save correlation vectors for figures.
#' 
#' OUTPUT
#' - `stats_plots.pdf`: PDF file containing histograms and scatter plots.
#' - `stats_plots.png`: PNG file converted from the PDF.
#' - `r_and_p.txt`: Text file containing correlation coefficients and p-values.
#' - `vectorsForFigs.RData`: RData file containing correlation vectors.
#' 
#' USAGE
#' Run the script from the command line with the required arguments:
#' ```
#' Rscript plot_meth_full.R <working_directory> <motif_length> <motif_symmetry> <score_cutoff>
#' ```
# ========================================

# ----- LOAD LIBRARY
args <- commandArgs(TRUE)
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(plyr, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(dplyr, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(scales, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(lazyeval, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(backports, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(corrplot, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(survival, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(Hmisc, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")
library(magrittr, lib.loc="/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/")

# ----- CONFIG

MINNUMPOINT=10

# ----- LOAD ARGUMENTS

motifLg=as.integer(args[2])
# print('Motif length:')
# print(motifLg)
SYM=args[3]
# print("Motif symetry:")
# print(SYM)
SC=as.integer(args[4]) # score cutoff to filter based on BBS
# print("Scores cut-off:")
# print(SC)
setwd(args[1])
vectors4Figs=paste(args[1], "vectorsForFigs.RData", sep="/")
# print(vectors4Figs)

# ----- FUNCTION
plot_and_test<-function(x,y,T,l,LI)
{
  #' This function performs a correlation test between two numeric vectors, `x` and `y`, 
  #' and generates a scatter plot with a regression line. It also updates a list with 
  #' the correlation coefficient and p-value.
  #'
  #' - x: A numeric vector representing the independent variable.
  #' - y: A numeric vector representing the dependent variable.
  #' - T: A character string used as the title prefix for the plot.
  #' - l: A character string used as the x-axis label for the plot.
  #' - LI: A list with elements `vR` and `vP` to store correlation coefficients 
  #'   and p-values, respectively.
  #'
  #' PROCESS
  #' The function filters out `NA` and `Inf` values from the input vectors `x` and `y`. 
  #' If the filtered vectors have more than a predefined minimum number of points 
  #' (`MINNUMPOINT`), it performs a Pearson correlation test and plots the data. 
  #' The plot includes a regression line and displays the correlation coefficient, 
  #' p-value, and the number of points used in the test. If the number of points is 
  #' insufficient, `NA` values are appended to the `LI` list.
  #'
  #' RETURN
  #' A list `LI` with updated elements:
  #'`vR`: A vector of correlation coefficients.
  #'`vP`: A vector of p-values.
  #' 
	x2=x[!is.infinite(y) & !is.na(x)]
    y2=y[!is.infinite(y) & !is.na(x)]
    if (length(x2)>MINNUMPOINT) {
    #test=cor.test(x2, y2, method="spearman", use="complete.obs")
    test=cor.test(x2, y2, method="pearson", use = "complete.obs")
    r=round(test$estimate, 2)
    P=round(test$p.value, 5)
	#P=round(-log10(test$p.value), 0)
    LI[["vR"]]=c(LI[["vR"]], r)
    LI[["vP"]]=c(LI[["vP"]], P)
    N=length(x2)
    plot(x2, y2, main=paste(T, "r:", r, "P:", P , "N:", N), 
       ylab="DAP/ampDAP signal ratio", xlab=l, cex=0.4, col="grey40", bg="grey40", lwd=1, pch=21, cex.lab=1.5, cex.axis=1.5, cex.main=1.7)
    abline(lm(y2 ~ x2), lwd=2, col="black")
    }
  else {
    LI[["vR"]]=c(LI[["vR"]], NA)
    LI[["vP"]]=c(LI[["vP"]], NA)
    }
  return(LI)
}

#-------------------------------
dodo<-function(N1, N2, SYM) {
  #' This function performs a series of operations based on the input parameters 
  #' and generates plots and statistical tests. It processes data from a given 
  #' dataset (`bs`) and applies different logic depending on whether the input 
  #' symmetry (`SYM`) is "yes" or "no".
  #'
  #' - N1 Character. Represents the nucleotide or motif to be analyzed.
  #' - N2 Character. Represents the complementary nucleotide or motif to be analyzed.
  #' - SYM Character. Indicates whether the analysis is symmetric ("yes") or not ("no").
  #'
  #' PROCESS
  #' - When `SYM` is "no", the function subsets the data for forward and reverse strands 
  #'   separately and performs operations on them.
  #' - When `SYM` is "yes", the function considers both strands together and applies 
  #'   symmetric logic for processing.
  #' - The function uses a loop to iterate over positions in the motif, calculates 
  #'   metrics, and generates plots using the `plot_and_test` function.
  #'
  #' RETURN
  #'  A list (`lili`) containing the results of the `plot_and_test` function 
  #' for each iteration.
  #'
  #' NOTE
  #' - The function assumes the existence of several global variables such as `motifLg`, 
  #'   `bs`, `S`, `xleg`, `lili`, `MINNUMPOINT`, and `cpt`.
  #' - The `plot_and_test` function is called within this function, and its behavior 
  #'   is critical to the output.

  if (SYM=="no") {size=motifLg-1}
  else if (SYM=="yes") {size=as.integer(motifLg/2)}
  
  for (i in 0:size) {
   #if (cpt%%9==0) {par(mfrow=c(3,3))}
   if (cpt==0) {par(mfrow=c(3,3))} # used for LFY and At1g19210
   if (SYM=="no") {
     sen<-subset(bs, bs[[S+i*2+1]]==N1 & bs$strand=="+")
     rev<-subset(bs, bs[[S+motifLg*2-i*2-1]]==N2 & bs$strand=="-")
     pm<-c(sen[[S+i*2+2]], rev[[S+motifLg*2-i*2]])
     r<-c(sen$r, rev$r)
	 lili<-plot_and_test(pm, r, paste(N1, i+1, sep=""), xleg, lili)
   }
   else if (SYM=="yes") {
    sel<-subset(bs, bs[[S+i*2+1]]==N1 | bs[[S+motifLg*2-i*2-1]]==N2) 
    pm<-pmax(sel[[S+i*2+2]], sel[[S+motifLg*2-i*2]], na.rm=TRUE)
    r<-sel$r
	lili<-plot_and_test(pm, r, paste(N1, i+1, N2, motifLg-i, sep=""), xleg, lili)
   }
   
   #lili<-plot_and_test(pm, r, paste(N1, i+1, N2, motifLg-i, sep=""), xleg, lili)
   #print(lili)
   if (length(pm)<=MINNUMPOINT) {cpt=cpt-1}
   cpt=cpt+1
  }
  return(lili)
}


# ----- MAIN 
lili=list("vR"=c(), "vP"=c()) #list initiation

# LOAD DATASET INPUT
d <-read.table("table_peak_bs_meth_Zhu.txt", header=T)

# head(d)
# summary(d$BSpfmScore) # TFBSs scores distribution

minDAP=min(d$covDAP[d$covDAP!=0])
# print(minDAP)
minAMP=min(d$covAMP[d$covAMP!=0])
# print(minAMP)

d2<- d %>% mutate(r=log10((covDAP+minDAP)/(covAMP+minAMP)), peakPropMeth=methSitesInPeak/(peakEnd-peakStart))

print("#########")
# head(d2)


# PLOTTING AND CALLS 
pdf("stats_plots.pdf", width=10, height=10)
hist(d$BSpfmScore, n=500)
par(mfrow=c(2,1))

## 1 - meth in peaks --------------------
#sel<-d2 %>% filter(EckerCovCGsInPeak/CGsInPeak>=0) %>% filter(peakPropMeth>=0)
sel<-d2 %>% filter(peakPropMeth>=0 & BSpfmScore>SC)
# print(nrow(sel))
lili<-plot_and_test(sel$peakPropMeth, sel$r, "Methylation in bound regions", "density of methylation (proportion of cytosine with a methylation probability >50%)", lili)

# writting list in a txt file -> Atlas heatmap representation
r_p <- paste(paste("r", gsub("[^0-9.-]", "", lili$vR), sep="\t"), paste("p", lili$vP, sep="\t"), sep="\n") #extract r (correlation)
fout="r_and_p.txt"
write.table(r_p, file=fout, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

bs <- d2 %>% filter(BSpfmScore>SC) %>% mutate(motifPropMeth=methSitesInMotif/motifLg) %>% mutate(sumPM = rowSums(.[grep("pm", names(.))], na.rm = T))
#write.table(bs, file=paste(args[3],"_bs_table.txt", sep=""), sep="\t", row.names = FALSE)

## 2 - meth in motif --------------------
selBS<-bs%>% filter(motifPropMeth>=0)
lili<-plot_and_test(selBS$methSitesInMotif, selBS$r, "Methylation in PFM best-scoring Binding Sites", "density of methylation (proportion of cytosine with a methylation probability >50%)", lili)

###### Per position regression, all pos both strand ---------
S=match("s1",names(d))-1
cpt=0
xleg="probability of cytosine methylation"

fPr<-sel$r
fPm<-sel$peakPropMeth
fBSr<-selBS$r
fBSm<-selBS$methSitesInMotif

lili<-dodo("C","G", SYM)
lili<-dodo("G","C", SYM)

fposR=lili[["vR"]][3:length(lili[["vR"]])]
fposP=lili[["vP"]][3:length(lili[["vP"]])]
save(fPr,fPm, fBSr, fBSm, fposR, fposP, file=vectors4Figs)
# print("corr vector lg is:")
# print(length(fposR))
# print(fposR)

## between positions meth prob correlation
# print("before filtering")
# head(d2)
tmp <- d2 %>% filter(strand=="+") %>% select(starts_with("pm"))
# print(head(tmp))
# print("end")
par(mfrow=c(1,1))

cor_5 <- rcorr(as.matrix(tmp), type="spearman")
print("-----")
print(cor_5)

M <- cor_5$r
print("-----")
print(M)

p_mat <- cor_5$P
print("-----")
print(p_mat)

p_mat[p_mat == 0] <- 1
print("-----")
print(p_mat)

p_mat2<-matrix(p.adjust(as.vector(p_mat), method="BH"),ncol=ncol(p_mat))
print("-----")
print(head(p_mat2))

# p_mat<-p_mat2
# corrplot(M, type = "upper", 
# 	          p.mat = p_mat, sig.level = 0.01,
# 		  main="")
# 		  #main="Methylation probability correlation between positions")

dev.off()

system("convert -density 144 -background white stats_plots.pdf stats_plots.png")



  