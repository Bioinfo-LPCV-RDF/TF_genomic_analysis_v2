#
# R Script: Coverage and Fold Change Analysis
#
# DESCRIPTION:
# This R script generates scatter plots to compare coverage ratios and fold changes between two samples.
# It visualizes the relationship between coverage in libraries (RiL) and peaks (RiP) and highlights
# statistically significant changes. 
#
# USAGE:
# The script is designed to be called from the command line with the following arguments:
#
# ```bash
# Rscript script_name.R <output_directory> <sample1_name> <sample2_name> <table_name> <type>
# ```
#
# ARGUMENTS:
# - `<output_directory>` : Directory where the plots will be saved.                                   
# - `<sample1_name>`     : Name of the first sample (column name in the table).                       
# - `<sample2_name>`     : Name of the second sample (column name in the table).                      
# - `<table_name>`       : Path to the input table file (TSV format).                                 
# - `<type>`             : Type of analysis (not used in the current script).                         
#
#
# OUTPUTS:
# The script generates several PNG files in the specified output directory:
# - `<sample1_name>_<sample2_name>_RiL.png`: Combined plot showing RiL coverage and fold change.
# - `<sample1_name>_<sample2_name>_RiP.png`: Combined plot showing RiP coverage and fold change.
# - `<sample1_name>_<sample2_name>_RiLandRiP.png`: Combined plot showing both RiL and RiP.
# - `<sample1_name>_<sample2_name>_RiPvsRiL.png`: Plot comparing RiP vs RiL.
# - `<sample1_name>_<sample2_name>_RiLstats.png`: Plot showing RiL with statistical significance.
# - `<sample1_name>_<sample2_name>_RiPstats.png`: Plot showing RiP with statistical significance.
#
# NOTES:
# - The script uses a logarithmic scale for both axes to better visualize the data distribution.
# - The script assumes that the input table contains columns named `<sample1_name>_RiL`, `<sample1_name>_RiP`, `<sample2_name>_RiL`, `<sample2_name>_RiP`, `logFC`, and `FDR`.


# LOAD LIBRARIES
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(ggplot2)
library(Cairo) 
library(cowplot)
library(scales)
# theme_set(theme_cowplot())

# LOAD ARGUMENTS
args=commandArgs(trailingOnly=TRUE)
out_dir=args[1]
name1=args[2]
name2=args[3]
tableName=args[4]
Type=args[5]

# CREATE OUTPUT DIR AND LOAD DATA
dir.create(file.path(out_dir),recursive=TRUE ,showWarnings = FALSE)
table <- read.table(tableName,sep="\t",header=TRUE,check.names=FALSE)
print(head(table))
table$name <- factor(table$name, levels = c(name1,name2,"both"))

# COMPUTATIONS
# compute coverage ratio
table$ratioL<-table[,paste(as.name(name1),"RiL",sep="_")] / (table[,paste(as.name(name2),"RiL",sep="_")])
table$ratioP<-table[,paste(as.name(name1),"RiP",sep="_")] / (table[,paste(as.name(name2),"RiP",sep="_")])

# PLOTTING
# RiL with MACs attributing colors
covRiL1<-ggplot(data=table, aes(y=table[,paste(as.name(name1),"RiL",sep="_")], x=table[,paste(as.name(name2),"RiL",sep="_")]))+
geom_point(alpha=0.4,aes(color=name)) +
scale_x_log10() + scale_y_log10() + theme_bw() + labs(color="peaks called \nby MACS3") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))

# RiL with CFC attributing colors
covRiL2<-ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiL",sep="_")],x=table[,paste(as.name(name2),"RiL",sep="_")]))+
geom_point(alpha=0.4,aes(color=log2(ratioL))) +
scale_x_log10() + scale_y_log10() + theme_bw() +
scale_color_gradientn(colours=rainbow(6),name=paste("CFC = \n", name1," / ",name2,sep=""),limits=c(-5,5),oob=squish) +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))

#  RiP with MACS attributing colors
covRiP1<-ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiP",sep="_")],x=table[,paste(as.name(name2),"RiP",sep="_")]))+
geom_point(alpha=0.4,aes(color=name)) +
scale_x_log10() + scale_y_log10() + theme_bw() + 
labs(color="peaks called \nby MACS3") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))

#  RiP with CFC attributing colors
covRiP2<-ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiP",sep="_")],x=table[,paste(as.name(name2),"RiP",sep="_")]))+
geom_point(alpha=0.4,aes(color= log2(ratioP))) +
scale_x_log10() + scale_y_log10() + theme_bw() +
scale_color_gradientn(colours=rainbow(6),name=paste("CFC = \n", name1," / ",name2,sep=""),limits=c(-5,5),oob=squish)+
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))


# COMPUTE 
table$FC=0
table$FC[table$logFC<=-1&table$FDR<=0.05]=-1
table$FC[table$logFC>=1&table$FDR<=0.05]=1
table$FC<-as.factor(table$FC)
table$logFC2=table$logFC
table$logFC2[table$FDR>0.05]=NA
print(summary(table))


# PLOTTING
covRiPstat<-ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiP",sep="_")],x=table[,paste(as.name(name2),"RiP",sep="_")]))+
geom_point(alpha=0.4,aes(color= logFC2)) +
scale_x_log10() + scale_y_log10() + theme_bw() +
scale_color_gradientn(colours=c("#39685D", "#39685D", "#F1D7AD" ,"#99545D","#99545D"),name="",limits=c(-5,5),oob=squish,na.value="#F1D7AD") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))

covRiLstat<-ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiL",sep="_")],x=table[,paste(as.name(name2),"RiL",sep="_")]))+
geom_point(alpha=0.4,aes(color= logFC2)) +
scale_x_log10() + scale_y_log10() + theme_bw() +
scale_color_gradientn(colours=c("#39685D", "#39685D", "#F1D7AD" ,"#99545D","#99545D"),name="",limits=c(-5,5),oob=squish,na.value="#F1D7AD") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,sep="")) +
ylab (paste("cov ",name1,sep=""))


# EXTRACT LEGEND
legend1 <-get_legend(covRiL1)
covRiL1 <- covRiL1 + theme(legend.position='none')
legend2 <-get_legend(covRiL2)
covRiL2 <- covRiL2 + theme(legend.position='none')
legend3 <-get_legend(covRiP1)
covRiP1 <- covRiP1 + theme(legend.position='none')
legend4 <-get_legend(covRiP2)
covRiP2 <- covRiP2 + theme(legend.position='none')

# EXPORTATION
Cairo(width = 1500, height = 500, file=paste(out_dir,"/",name1,"_",name2,"_RiL.png",sep=""), type="png")
ggdraw() +
draw_plot(covRiL1, 0.0, 0.0, 0.5, 1) +
draw_plot(covRiL2, 0.5, 0, 0.5, 1) +

draw_plot(legend1, 0.43, -0.02, 0.2, 0.5)+
draw_plot(legend2, 0.90, -0.02, 0.18, 0.5)+

draw_plot_label(c("Reads mapped in library", "Reads mapped in library"), c(0.02, 0.52), c(0.995, 0.995), size = 15)
dev.off()

Cairo(width = 1500, height = 500, file=paste(out_dir,"/",name1,"_",name2,"_RiP.png",sep=""), type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
ggdraw() +
draw_plot(covRiP1, 0.0, 0.0, 0.5, 1) +
draw_plot(covRiP2, 0.5, 0, 0.5, 1) +

draw_plot(legend3, 0.43, -0.02, 0.2, 0.5)+
draw_plot(legend4, 0.90, -0.02, 0.18, 0.5)+

draw_plot_label(c("Reads mapped in peaks", "Reads mapped in peaks"), c(0.02, 0.52), c(0.995, 0.995), size = 15)
dev.off()

Cairo(width = 1500, height = 1000, file=paste(out_dir,"/",name1,"_",name2,"_RiLandRiP.png",sep=""), type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
ggdraw() +
draw_plot(covRiP1, 0.0, 0.5, 0.5, 0.5) +
draw_plot(covRiP2, 0.5, 0.5, 0.5, 0.5) +
draw_plot(covRiL1, 0.0, 0.0, 0.5, 0.5) +
draw_plot(covRiL2, 0.5, 0, 0.5, 0.5) +

draw_plot(legend1, 0.43, 0.48, 0.2, 0.25)+
draw_plot(legend4, 0.90, 0.48, 0.18, 0.25)+
draw_plot(legend3, 0.43, -0.02, 0.2, 0.25)+
draw_plot(legend2, 0.90, -0.02, 0.18, 0.25)+

draw_plot_label(c("Reads mapped in peaks", "Reads mapped in peaks"), c(0.02, 0.52), c(0.995, 0.995), size = 15)+
draw_plot_label(c("Reads mapped in library", "Reads mapped in library"), c(0.02, 0.52), c(0.495, 0.495), size = 15)
dev.off()

covPvsL1<-ggplot(data=table,aes(y=table[,paste(name1,"RiP",sep="_")],x=table[,paste(name1,"RiL",sep="_")]))+
geom_point(alpha=0.4,aes(color=name)) +
scale_x_log10() +
scale_y_log10() +
labs(color="peaks called \nby MACS2") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name1,"\nnormalized with peaks mapped in library",sep="")) +
ylab (paste("cov ",name1,"\nnormalized with peaks mapped in peaks called",sep=""))

covPvsL2<-ggplot(data=table,aes(y=table[,paste(name2,"RiP",sep="_")],x=table[,paste(name2,"RiL",sep="_")]))+
geom_point(alpha=0.4,aes(color=name)) +
scale_x_log10() +
scale_y_log10() +
labs(color="peaks called \nby MACS2") +
theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
geom_abline(intercept = 0) +
xlab (paste("cov ",name2,"\nnormalized with peaks mapped in library",sep="")) +
ylab (paste("cov ",name2,"\nnormalized with peaks mapped in peaks called",sep=""))


legend5 <-get_legend(covPvsL1)
covPvsL1 <- covPvsL1 + theme(legend.position='none')

legend6 <-get_legend(covPvsL2)
covPvsL2 <- covPvsL2 + theme(legend.position='none')

Cairo(width = 1500, height = 500, file=paste(out_dir,"/",name1,"_",name2,"_RiPvsRiL.png",sep=""), type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
# plot_grid(covPvsL1,covPvsL2, nrow=1, ncol=2, labels=c('A','B'))
ggdraw() +
draw_plot(covPvsL1, 0.0, 0.0, 0.5, 1) +
draw_plot(covPvsL2, 0.5, 0, 0.5, 1) +

draw_plot(legend5, 0.43, -0.02, 0.2, 0.5)+
draw_plot(legend6, 0.90, -0.02, 0.18, 0.5)+

draw_plot_label(c(paste(name1,": Reads mapped in peaks / Reads mapped in library"), paste(name2,": Reads mapped in peaks / Reads mapped in library")), c(-0.03, 0.47), c(0.995, 0.995), size = 15)
dev.off()

Cairo(width = 750, height = 750, file=paste(out_dir,"/",name1,"_",name2,"_RiLstats.png",sep=""), type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
covRiLstat
dev.off()
Cairo(width = 750, height = 750, file=paste(out_dir,"/",name1,"_",name2,"_RiPstats.png",sep=""), type="png", pointsize=12, bg = "white", canvas = "white", units = "px", dpi = "auto")
covRiPstat
dev.off()

quit(save = "no", status = 0, runLast = FALSE)