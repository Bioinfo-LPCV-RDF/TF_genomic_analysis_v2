#!/usr/bin/env Rscript

# Script: figs_meth_violin.R
# ===============================
# DESCRIPTION:
#   This script generates a multi-panel figure combining:
#     - PFM logo 
#     - Heatmaps of correlation data
#     - Scatter plots and violin plots for methylation signal data
# 
# RELATED TO cooking_methylation
# 
# USAGE:
#   Rscript methylation_tfbs_plot.R <working_directory> <pfm_file> <half_logo>
#
# ARGUMENTS:
#   working_directory : Path to the directory where input files are located and output will be saved
#   pfm_file          : Path to the monomer PFM file (tab-delimited, 4 columns: A C G T)
#   half_logo         : "yes" or "no", whether to use only the first half of the PFM for the logo
#
# Notes:
#   - Input files like 'vectorsForFigs.RData' must be present in the working directory.
#   - Heatmaps, scatter plots, violin plots, and sequence logos are integrated into a single figure.
#   - Output files: Figure.pdf and Figure.png in the working directory.

# LOAD ARGUMENTS 
args <- commandArgs(TRUE)

pfm=args[2]
sym=args[3]

# LOAD LIBRARIES
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(dplyr)
library(corrplot)
library(cowplot)
library(ggseqlogo)
library(ggplot2)
library(multcompView)

# - - - - - Preparation of methylation LOGOs - - - - - - - - -
pfm2logo<-function(PFM, SYM){
  system(paste('tail -n +2', PFM, "> pfm.txt", sep="\t"))
  pfm<-read.table("pfm.txt", header=T)
  if (SYM=="yes"){pfm<-head(pfm, ceiling(nrow(pfm)/2))} #if motif is symetric, cut half the logo
  logo=ggseqlogo(t(pfm))
  return(logo)
}

logo=pfm2logo(pfm, sym) # call to get the logo from PFM 

# - - - - -  Methylation - - - - -

# - - - - HEATMAP - - - - 
myHeatMap<-function(R, P) {
  L=length(R)/2
  Var1=c(c(c(1:L),c(1:L)))
  Var2=c(rep("5'", L), rep("3'", L))
  Pa<-p.adjust(P, method="BH")
  print(Pa)
  
  R=as.vector(R)
  R[Pa>0.05]<-0
  
  tmp=data.frame(cbind(as.character(Var1), Var2, as.numeric(R), as.numeric(Pa)))
		
  hm<-ggplot(data = tmp, aes(x=Var1, y=Var2, fill=R)) +  # size=log10(Pa)
	geom_tile(colour="grey") +
	labs(x=NULL, y=NULL) +
	scale_fill_gradient2(low = "#6D9EC1", high = "#E46726", mid = "white", 
		midpoint = 0, limit = c(-1,1), 
		name="r", na.value="grey80") +
		scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))+
		theme_void()

  hm <- hm + theme(legend.position="bottom",
		  legend.title=element_text(vjust=0.8),
		  axis.text.x=element_blank(),
		  panel.grid.major = element_blank(), 
		  panel.grid.minor = element_blank(), 
		  panel.background = element_blank(),
		  axis.ticks = element_blank(),
		  panel.border = element_blank())
  hm <- hm + theme(plot.margin=unit(c(t=0,r=0,b=0,l=0.),"cm"), axis.text.y = 		element_text(vjust = .5, hjust = 0))
  return(hm)
}



  
# - - - - SCATTER PLOTS - - - - - 
myPlot<-function(X, Y, XLAB, MINY, MAXY, XDOT, TITLE) {
  d<-data.frame(X, Y)
  t=cor.test(X, Y, method="pearson", use = "complete.obs")
  r=round(t$estimate, 2)
  #P=round(t$p.value, 10)
  P=round(-log10(t$p.value), 0)
  if (length(XDOT)>1) {
    ran<-sample(seq(from=-0.2, to=0.2, by=0.01), replace=TRUE, size=length(Y))
    X2=X+ran
  }
  else {X2=X}
  d<-data.frame(X2, Y)
  p<-ggplot(d, aes(x=X2, y=Y)) +
       geom_point(colour="grey40", size=0.1) +
       stat_smooth(method=lm, colour="black") +
	   #ylim(MINY, MAXY)+
       xlab(XLAB) + 
	   #ylab("DAP/ampDAP signal ratio") +
	   ylab(bquote(~log[10]~ '(DAP/ampDAP signal)')) +
	   theme(axis.title=element_text(size=10), 
	   plot.margin = margin(t=7, r=10, b=7, l=7))+
       annotate('text', max(X), max(Y),
			label=paste("italic(r)==", r, '~-log[10](', "~italic(P))==", P), parse=TRUE, hjust=1, size=3)
	   
       if (length(XDOT)>1) {
	    p<-p + scale_x_discrete(limits=XDOT)
       }
	   if (TITLE!="NO") {
	     p<-p + ggtitle(TITLE) + 
		 theme(plot.title = element_text(hjust = 0.5, size=10))}
  return(p)
}

# - - - - VIOLIN PLOTS - - - - - 
myViolin<-function(X, Y, XLAB, MINY, MAXY, XDOT, TITLE) {
  d<-data.frame(X, Y)
  t=cor.test(X, Y, method="pearson", use = "complete.obs")
  r=round(t$estimate, 2)
  P=round(-log10(t$p.value), 0)
  cf <- coef(lm(Y~X, data=d))

  d<-data.frame(X, Y)
  d$X <- as.factor(X)
  
  # remove categories with less than 3 observations
  counts <- table(d$X)
  d_subset <- d[d$X %in% names(counts[counts >= 3]), ]
  X <- d_subset$X
  Y <- d_subset$Y
  d<-data.frame(X, Y)
  d$X <- as.factor(X)
  
  aov_model <- aov(Y ~ X, data = d)
  summary(aov_model)
  tukey <- TukeyHSD(aov_model, ordered = TRUE)
  print(tukey)
  cld <- multcompLetters4(aov_model, tukey)
  print(cld)
  
  # table with factors and 3rd quantile
  Tk <- group_by(d, X) %>%
    summarise(mean=mean(Y), quant = quantile(Y, probs = 0.9)) %>%
    arrange(desc(mean))
	
  # extracting the compact letter display and adding to the Tk table
  cld <- as.data.frame.list(cld$X)
  Tk$cld <- cld$Letters
  print(Tk)
  
  # Extract pairwise comparisons and format as a string
  #tukey_comparisons <- capture.output(print(tukey_test))
  #tukey_comparisons <- paste(tukey_comparisons, collapse = "\n")
  

  p<-ggplot(d, aes(x=X, y=Y)) +
       geom_violin(fill="#ffdce2", colour="black", size=0.8) +
	   
	   #stat_summary(fun.data="mean_sdl", mult=1, 
       #          geom="pointrange", color="black" )+
	   stat_summary(aes(group=1),fun=mean, geom="point", color="black", size=1.2) + 
	   #geom_abline(slope=cf[2], intercept=cf[1], lwd=.4) +
       #stat_smooth(method=lm, colour="black") +
	   #ylim(MINY, MAXY)+
       xlab(XLAB) +  ylab(bquote(~log[10]~ '(DAP/ampDAP signal)')) +
	   theme(axis.title=element_text(size=10)) +
	   geom_text(data = Tk, aes(x = X, y = quant, label = cld ), size = 3, vjust=-1, hjust =-1)
	   
	   if (TITLE!="NO") {
	     p<-p + ggtitle(TITLE) + 
		 theme(plot.title = element_text(hjust = 0.5, size=10))}
  return(p)
}



# - - - - Plot the figure - - - -
setwd(args[1])
load("vectorsForFigs.RData") # assumes vectors fposR, fposP, fPm, fPr, fBSm, fBSr exist
hm=myHeatMap(fposR, fposP)
m=0.1
minY=fPr-m
maxY=fPr+m
fp<-myPlot(fPm, fPr, "Methylation density", minY, maxY, c(), "")

minY=fBSr-m
maxY=fBSr+m
fm<-myViolin(fBSm, fBSr, "#methylated cytosines in TF binding motif", minY, maxY, c(0,1,2,3,4,5,6,7,8), "NO")


pdf("Figure.pdf", width=4, height=7)
ggdraw() +
	  draw_plot(fp, 0, .66, 1, 0.33) +
	  draw_plot(fm, 0, .33, 1, .33) +
	  draw_plot(logo, .0, .17, 1, .17) +
	  draw_plot(hm, x = 0.125, y = .065, width =.81 , height = 0.13)
dev.off()

system("convert -density 144 -background white Figure.pdf Figure.png")

