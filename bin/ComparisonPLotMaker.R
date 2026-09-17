# FUNCTION
#  Comparison Plot Maker
# ======================
# 
# RELATED TO pairwize_comparison
#  
# DESCRIPTION
#  This script generates various comparison plots for genomic analysis using ggplot2. 
#  It processes input data, computes metrics, and creates visualizations to compare two datasets.
# 
# USAGE Rscript ComparisonPlotMaker.R -t <tableName> -od <outdir> -n1 <name1> -n2 <name2> 
#        -c1 <color1> -c2 <color2> -cn <colorneutral> [-th <thresholdFC>]
# 
# ARGUMENTS
# -t, --tableName Path to the input table file (tab-delimited, with headers).
# -od, --outdir Output directory where plots will be saved.
# -n1, --name1 Name of the first dataset for comparison.
# -n2, --name2 Name of the second dataset for comparison.
# -c1, --color1 Color for the first dataset in plots.
# -c2, --color2 Color for the second dataset in plots.
# -cn, --colorneutral Neutral color for overlapping data points.
# -th, --thresholdFC (Optional) Threshold for fold change (default: 1).
# 
# WORKFLOW
# The script performs the following tasks:
# - Reads the input table and processes column names to identify data for plotting.
# - Generates scatter plots comparing coverage between two datasets.
# - Computes Coverage Fold Change (CFC) and visualizes it in various contexts.
# - Uses EdgeR statistics (if available) to identify peaks differentially covered (PDC).
# - Saves plots as PNG files in the specified output directory.
# 
# NOTE
# - The script assumes specific column naming conventions in the input table.
# - Requires the following R packages: ggplot2, Cairo, cowplot, scales, plotROC, RColorBrewer, argparse.
# - Ensure the input table contains the necessary columns for the selected suffixes (e.g., "RiLuc", "RiL", "RiP").
# 

.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(ggplot2)
library(Cairo) 
library(cowplot)
library(scales)
library(plotROC)
library(RColorBrewer)
library("argparse")

parser <- ArgumentParser()

parser$add_argument("-t", "--tableName")
parser$add_argument("-od", "--outdir")
parser$add_argument("-n1", "--name1")
parser$add_argument("-n2", "--name2")
parser$add_argument('-c1', "--color1")
parser$add_argument('-c2', "--color2")
parser$add_argument('-cn', "--colorneutral")
parser$add_argument('-th', "--thresholdFC", default=1, type="double")

args <- parser$parse_args()
tableName <- args$tableName
out_dir <- args$outdir
name1 <- args$name1
name2 <- args$name2
color1 <- args$color1
color2 <- args$color2
colorneut <- args$colorneutral
thresholdFC <- args$thresholdFC

#Creating directory if need be
dir.create(file.path(out_dir),recursive=TRUE ,showWarnings = FALSE)
table <- read.table(tableName,sep="\t",header=TRUE,check.names=FALSE)
# getting the colnames to identify what to plot
tablecolnames<-colnames(table)

table$name <-factor(table$name, levels=c(name1,name2,"both"))

colfunc <- colorRampPalette(c(color1, colorneut, color2))

plots_list<-list()
plotsname_v<-c()
plotslabel_v<-c()
nbplots=1

prefixplotsname=paste(out_dir,"/",name1,"_",name2,sep="")
summary(table)
vec_suffix=c("RiLuc","RiL","RiP")
# vec_suffix=c("RiP")
# looping on all suffixes possible
for (suffix in vec_suffix){
	# checking if colnames exists to prepare plots
	if((paste(name1,suffix,sep="_") %in% tablecolnames)&(paste(name2,suffix,sep="_") %in% tablecolnames)){
		print(paste("preparing plots for",suffix))
		# plot with MACS3 determined peaks

		# summary(table[,paste(as.name(name1),suffix,sep="_")])
		print(nbplots)
		print(paste(as.name(name1),suffix,sep="_"))
		tmp1 <- ggplot(data=table) +
		geom_point(alpha=0.4,aes(color=name,y=table[,paste(as.name(name1),suffix,sep="_")], x=table[,paste(as.name(name2),suffix,sep="_")])) +
		scale_x_log10() + scale_y_log10() + theme_bw() +
		labs(color="peaks called \nby MACS3") +
		theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
		scale_color_manual(labels=c(name1,name2,"both"), values = c(color1, color2, colorneut)) +
		geom_abline(intercept = 0) +
		guides(colour = guide_legend(override.aes = list(size=2.5))) +
		xlab (paste("cov",name2,suffix,sep=" ")) +
		ylab (paste("cov",name1,suffix,sep=""))
		ggsave2(filename=paste(prefixplotsname,"_MACS3_",suffix,".png",sep=""),plot=tmp1,width=7,height=5,units="in")

		# computing CFC (Coverage Fold Change)
		table[,paste("CFC",suffix,sep="_")] <- table[,paste(as.name(name1),suffix,sep="_")] / (table[,paste(as.name(name2),suffix,sep="_")])
		# plots with CFC determined peaks
		tmp2 <- ggplot(data=table,aes(y=table[,paste(as.name(name1),suffix,sep="_")],x=table[,paste(as.name(name2),suffix,sep="_")]))+
		geom_point(alpha=0.4,aes(color=log2(table[,paste("CFC",suffix,sep="_")]))) +
		scale_x_log10() + scale_y_log10() + theme_bw() +
		scale_color_gradientn(colours=colfunc(6),name=paste("CFC = \n", name1," / ",name2,sep=""),limits=c(-5,5),oob=squish) +
		theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
		geom_abline(intercept = 0) +
		xlab (paste("cov",name2,suffix,sep=" ")) +
		ylab (paste("cov",name1,suffix,sep=""))
		ggsave2(filename=paste(prefixplotsname,"_CFC_",suffix,".png",sep=""),plot=tmp2,width=7,height=5,units="in")
		
		#log2CFC in function of mean coverage
		tmp2b <- ggplot(data=table,aes(x=(table[,paste(as.name(name1),suffix,sep="_")]+table[,paste(as.name(name2),suffix,sep="_")])/2,y=log2(table[,paste("CFC",suffix,sep="_")]))) +
		geom_point(alpha=0.4,aes(color=log2(table[,paste("CFC",suffix,sep="_")]))) +
		scale_x_log10() + theme_bw() +
		scale_color_gradientn(colours=colfunc(6),name=paste("CFC = \n", name1," / ",name2,sep=""),limits=c(-5,5),oob=squish) +
		theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
		geom_hline(yintercept = 0) +
		ylim(-10,10) +
		xlab (paste("Mean cov",suffix,sep=" ")) +
		ylab (paste("log2 CFC",name1,"/",name2,sep=" "))
		ggsave2(filename=paste(prefixplotsname,"_log2CFC_vsMeanCov_",suffix,".png",sep=""),plot=tmp2b,width=7,height=5,units="in")
		# log2CFC in function of coverage of name1
		tmp2c <- ggplot(data=table,aes(x=table[,paste(as.name(name1),suffix,sep="_")],y=log2(table[,paste("CFC",suffix,sep="_")]))) +
		geom_point(alpha=0.4,aes(color=log2(table[,paste("CFC",suffix,sep="_")]))) +
		scale_x_log10() + theme_bw() +
		scale_color_gradientn(colours=colfunc(6),name=paste("CFC = \n", name1," / ",name2,sep=""),limits=c(-5,5),oob=squish) +
		theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
		geom_hline(yintercept = 0) +
		ylim(-10,10) +
		xlab (paste("cov",name1,suffix,sep=" ")) +
		ylab (paste("log2 CFC",name1,"/",name2,sep=" "))
		ggsave2(filename=paste(prefixplotsname,"_log2CFC_vsCov_",name1,"_",suffix,".png",sep=""),plot=tmp2c,width=7,height=5,units="in")
		# log2CFC in function of coverage of name2
		tmp2d <- ggplot(data=table,aes(x=table[,paste(as.name(name2),suffix,sep="_")],y=log2(table[,paste("CFC",suffix,sep="_")]))) +
		geom_point(alpha=0.4,aes(color=log2(table[,paste("CFC",suffix,sep="_")]))) +
		scale_x_log10() + theme_bw() +
		scale_color_gradientn(colours=colfunc(6),name=paste("CFC = \n", name1," / ",name2,sep=""),limits=c(-5,5),oob=squish) +
		theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
		geom_hline(yintercept = 0) +
		ylim(-10,10) +
		xlab (paste("cov",name2,suffix,sep=" ")) +
		ylab (paste("log2 CFC",name1,"/",name2,sep=" "))
		ggsave2(filename=paste(prefixplotsname,"_log2CFC_vsCov_",name2,"_",suffix,".png",sep=""),plot=tmp2d,width=7,height=5,units="in")

		
		# Using EdgeR stats to detect PDC (peaks differentially covered; workd like DEG for RNA-seq analysis)
		if( (suffix == "RiP")&("FDR" %in% tablecolnames) ){
			print("using stats")
			# peaks with high FC (in a way or another), but without a good FDR (meaning their coverages varies too much in between replicates of same sample) are tagged "NA"
			table$FC="both"
			table$FC[table$logFC<(-thresholdFC)&table$FDR<=0.05]=name1
			table$FC[table$logFC>thresholdFC&table$FDR<=0.05]=name2
			# table$FC[table$logFC<1&table$logFC>-1]="both"
			table$FC<-factor(table$FC,levels=c("both",name1,name2))
			print(summary(table))
			# plots with EdgeR determined peaks
			tmp3 <- ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiP",sep="_")],x=table[,paste(as.name(name2),"RiP",sep="_")]))+
			geom_point(alpha=0.4,aes(color= FC)) +
			scale_x_log10() + scale_y_log10() + theme_bw() +
			theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
			geom_abline(intercept = 0) +
			guides(colour = guide_legend(override.aes = list(size=2.5))) +
			scale_color_manual(labels=c("both", name1,name2,'NA'), values = c(colorneut, color1, color2,"#000000")) +
			xlab (paste("cov",name2,suffix,sep=" ")) +
			ylab (paste("cov",name1,suffix,sep=""))
			ggsave2(filename=paste(prefixplotsname,"_EdgeR_",suffix,".png",sep=""),plot=tmp3,width=7,height=5,units="in")

			tmp3 <- ggplot(data=table,aes(y=table[,paste(as.name(name1),"RiP",sep="_")],x=table[,paste(as.name(name2),"RiP",sep="_")]))+
			geom_point(alpha=0.4,aes(color= FDR)) +
			scale_x_log10() + scale_y_log10() + theme_bw() +
			theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
			geom_abline(intercept = 0) +
			# guides(colour = guide_legend(override.aes = list(size=2.5))) +
			scale_color_gradientn(colours=rev(rainbow(6)),name=paste("FDR : \n", name1," vs ",name2,sep=""),limits=c(0,0.25),oob=squish) +
			xlab (paste("cov",name2,suffix,sep=" ")) +
			ylab (paste("cov",name1,suffix,sep=""))
			ggsave2(filename=paste(prefixplotsname,"_FDR_",suffix,".png",sep=""),plot=tmp3,width=7,height=5,units="in")

			# log2CFC in function of mean coverage with EdgeR determined peaks (colored by FDR)
			tmp3b <- ggplot(data=table,aes(x=(table[,paste(as.name(name1),suffix,sep="_")]+table[,paste(as.name(name2),suffix,sep="_")])/2,y=log2(table[,paste("CFC",suffix,sep="_")]))) +
			geom_point(alpha=0.4,aes(color= FDR)) +
			scale_x_log10() + theme_bw() +
			scale_color_gradientn(colours=rev(rainbow(6)),name=paste("FDR : \n", name1," vs ",name2,sep=""),limits=c(0,0.25),oob=squish) +
			theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
			ylim(-10,10) +
			geom_hline(yintercept = 0) +
			xlab (paste("Mean cov",suffix,sep=" ")) +
			ylab (paste("log2 CFC",name1,"/",name2,sep=" "))
			ggsave2(filename=paste(prefixplotsname,"_log2CFC_vsMeanCov_EdgeR_",suffix,".png",sep=""),plot=tmp3b,width=7,height=5,units="in")

			# log2CFC in function of coverage of name1 with EdgeR determined peaks (colored by FDR)
			tmp3c <- ggplot(data=table,aes(x=table[,paste(as.name(name1),suffix,sep="_")],y=log2(table[,paste("CFC",suffix,sep="_")]))) +
			geom_point(alpha=0.4,aes(color= FDR)) +
			scale_x_log10() + theme_bw() +
			scale_color_gradientn(colours=rev(rainbow(6)),name=paste("FDR : \n", name1," vs ",name2,sep=""),limits=c(0,0.25),oob=squish) +
			theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
			ylim(-10,10) +
			geom_hline(yintercept = 0) +
			xlab (paste("cov",name1,suffix,sep=" ")) +
			ylab (paste("log2 CFC",name1,"/",name2,sep=" "))
			ggsave2(filename=paste(prefixplotsname,"_log2CFC_vsCov_EdgeR_",name1,"_",suffix,".png",sep=""),plot=tmp3c,width=7,height=5,units="in")

			# log2CFC in function of coverage of name2 with EdgeR determined peaks (colored by FDR)
			tmp3d <- ggplot(data=table,aes(x=table[,paste(as.name(name2),suffix,sep="_")],y=log2(table[,paste("CFC",suffix,sep="_")]))) +
			geom_point(alpha=0.4,aes(color= FDR)) +
			scale_x_log10() + theme_bw() +
			scale_color_gradientn(colours=rev(rainbow(6)),name=paste("FDR : \n", name1," vs ",name2,sep=""),limits=c(0,0.25),oob=squish) +
			theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
			ylim(-10,10) +
			geom_hline(yintercept = 0) +
			xlab (paste("cov",name2,suffix,sep=" ")) +
			ylab (paste("log2 CFC",name1,"/",name2,sep=" "))
			ggsave2(filename=paste(prefixplotsname,"_log2CFC_vsCov_EdgeR_",name2,"_",suffix,".png",sep=""),plot=tmp3d,width=7,height=5,units="in")
		}
	}
}
# for(name in c(name1,name2)){
# 	if((paste(name,'RiLuc',sep="_") %in% tablecolnames)&(paste(name,'RiP',sep="_") %in% tablecolnames)){
# 		tmp4 <- ggplot(data=table, aes(y=table[,paste(as.name(name),'RiLuc',sep="_")], x=table[,paste(as.name(name),'RiP',sep="_")])) +
# 		geom_point(alpha=0.4,aes(color=name)) +
# 		scale_x_log10() + scale_y_log10() + theme_bw() +
# 		labs(color="peaks called \nby MACS3") +
# 		theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
# 		geom_abline(intercept = 0) +
# 		guides(colour = guide_legend(override.aes = list(size=2.5))) +
# 		xlab (paste("cov ",name,' RiP',sep="")) +
# 		ylab (paste("cov ",name,' RiL',sep=""))
# 		ggsave2(filename=paste(prefixplotsname,"RiPvsRiL",name,"test.png",sep="_"),plot=tmp4,width=7,height=5,units="in")
# 	}
# }
### Tests de détermination de groupes basé sur les ROCAUC
# nbrep=length(table[,paste(as.name(name1),"RiP",sep="_")])
# df <- data.frame(d=c(rep(1,nbrep),rep(0,nbrep)),M1=c(table[,paste(as.name(name1),"RiP",sep="_")],table[,paste(as.name(name2),"RiP",sep="_")]))
# basicplot <- ggplot(df, aes(d = d, m = M1)) + geom_roc(n.cuts=0)
# basicplot <- basicplot + annotate("text", x = .75, y = .25, 
#            label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2)))
# ggsave2(filename=paste(prefixplotsname,name,"ROC","test.png",sep="_"),plot=basicplot,width=7,height=5,units="in")
