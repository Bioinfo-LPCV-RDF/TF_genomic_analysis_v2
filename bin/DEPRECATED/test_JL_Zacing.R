#!/usr/bin/env Rscript

#Rscript Zacing.R SEP3AGspe/SEP3AGspe_spacing_pos.tsv 0.4 yes 0.3 0.4 0.5
#Rscript Zacing.R ARF5_spacing_pos.tsv no -4 -6 -7

# - - - set arguments - - -
#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
	make_option(c("-f", "--file"), type="character", default=NULL, 
		help="Spacing recap file", metavar="character"),
	make_option(c("-o", "--outdir"), type="character", default=NULL, 
		help="Spacing recap file", metavar="character"),
    make_option(c("-m", "--mat_type"), type="character", default="ASYMMETRIC", 
		help="matrix type [default %default]", metavar="character"),
	make_option(c("-t", "--thresholds"), type="character", default="[-8,-9,-10]", 
		help="thresholds for matrix, separate them with a single comma [default %default]"),
	make_option(c("-s", "--FCthresholds"), type="double", default=0.6, 
		help="thresholds for FoldChange [default %default]")
	); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (opt$mat_type!="ASYMMETRIC" && opt$mat_type!="SYMMETRIC"){quit(save="no", status=1)}

opt$thresholds=substr(opt$thresholds,2,nchar(opt$thresholds)-1)
opt$thresholds=strsplit(opt$thresholds, ",")

meth="median" #the code accepts mean and uniform but mean is bad. uniform gives similar resultat as median.
setwd(opt$outdir)

if (file.exists("Zscore_summary.txt")){
	system("rm Zscore_summary.txt")
}


#library(plyr)
#library(scales)
#library(lazyeval)
#library(backports)

#options(bitmapType='cairo')
library(dplyr,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
# library(tidyverse)
library(cowplot,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
theme_set(theme_cowplot())
library(ggrepel,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)





# - - - define functions - - -
foldchange<-function(BP, CONF, DF, S) {
	population=nrow(DF)
	if (CONF=="AR"){
		DF<-DF
	} else {
		DF<-DF %>% filter(Conf==CONF)
	}
	x<-nrow(DF%>%filter(Space==BP))
	Scount<-table(DF$Space) # this is to compute mean and median
	
	if (meth=="uniform") {
		
		u=population*(1/(length(S)*3)) # same as population * (max - min)/2
	}
	else if (meth=="mean") {u=mean(Scount)}
	else if (meth=="median") {u=median(Scount)}
	
	FC=x/u
	#print("fold change")
	#print(FC)
	return(FC)	
}

Zscore_calc<-function(BP, CONF, DF, S) {
#formula adpated from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142402/
# but finaly not really adapted from the article above...;
	population=nrow(DF)
	if (CONF=="AR"){
		DF<-DF
	} else {
		DF<-DF %>% filter(Conf==CONF)
	}
	
	x<-nrow(DF%>%filter(Space==BP))
	Scount<-table(DF$Space)
	
	if (meth=="uniform") {
		
		u=population*(1/(length(S)*3)) # same as population * (max - min)/2
		#print("sigma")
		P=1/(length(S)*3)
		sigma=sqrt(population*P*(1-P))
		#print(sigma)
		#OR (not exactly the same values obtained...)
		#variance=((max(S)-min(S))^2)/12
		#sigma=sqrt(variance)
		#print(sigma)
	}
	# compute the real u from mean or median ♫
	else if (meth=="mean") {
		u=mean(Scount)
		sigma=sd(Scount)
	}
	else if (meth=="median") {
		u=median(Scount)
		sigma=mad(Scount) # mad computes the median absolute deviation
	}
	
	Z=(x-u)/sigma	
	return(Z)
}

Zscore2Pvalue<-function(myZ){
	p=(1/sqrt(2*3.141593)*exp((-myZ^2)/2)) # formula taken from here: https://goodcalculators.com/p-value-calculator/
    return(p)
}

convert.z.score<-function(z, one.sided=NULL) {
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


myplot<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC) {
	D<-data.frame(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC)
	
	
	D$tmp<-NA
	D$tmp[abs(log2(D$FC))>opt$FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>opt$FCthreshold & D$ADJUSTED_PVAL<0.05]
	# with the median method I could aslo filter by Zscore> 3.5 as recommended by Iglewicz, B. and Hoaglin, D.C. (2010)
	
	res<-D$tmp[!is.na(D$tmp)]
	write(res,file="Zscore_summary.txt",append=TRUE)
	
	p<-ggplot(D, aes(x=SPA, y=ZCO, label=tmp)) +
       geom_point(colour="grey40", size=2.5) + geom_line() +
	   geom_text_repel() +
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}

volcanoPlot<-function(SPA,FC,PVAL, ADJUSTED_PVAL) {
	
	D<-data.frame(SPA, FC, PVAL, ADJUSTED_PVAL)
	
	D$tmp<-NA
	D$tmp[abs(log2(D$FC))>opt$FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>opt$FCthreshold & D$ADJUSTED_PVAL<0.05]
	myCol="darkblue"
	
	
	p<-ggplot(D,aes(x=log2(FC), y=-log10(ADJUSTED_PVAL), label=tmp)) + 
		geom_point() + 
		geom_text_repel() +
		xlab("log2(FC)") + ylab("-log10(P)") +
		geom_vline(xintercept=c(-opt$FCthreshold, opt$FCthreshold), col=myCol, linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), col=myCol, linetype="dashed")
	return(p)
}


main<-function(DATA, SCO, F) {
	
	write("#- - - - - - - - - - - - - - - - -",file="Zscore_summary.txt",append=TRUE)
	txt=paste(F, SCO, sep="_")
	write(txt,file="Zscore_summary.txt",append=TRUE)
	filtDF<-DATA%>%filter(Score1>SCO & Score2>SCO) 
# 	print(head(filtDF))
# 	print(summary(filtDF))
	
	minSpa=min(as.numeric(filtDF$Space))
	maxSpa=max(as.numeric(filtDF$Space))
	print(minSpa)
	print(maxSpa)
	#quit()
	Spacings<-seq(minSpa,maxSpa)
	outfig=paste("Zscore_", F, ".pdf", sep="")
	outtable=paste("Zscore_stats_", F, ".tsv", sep="")
	
	if (opt$mat_type=="ASYMMETRIC") {
		
		FC_ER<-sapply(Spacings, foldchange, DF=filtDF, CONF="ER", S=Spacings)
		FC_IR<-sapply(Spacings, foldchange, DF=filtDF, CONF="IR", S=Spacings)
		FC_DR<-sapply(Spacings, foldchange, DF=filtDF, CONF="DR", S=Spacings)
	
		Z_ER<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="ER", S=Spacings)
		Z_IR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="IR", S=Spacings)
		Z_DR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="DR", S=Spacings)
	
		P_ER<-sapply(Z_ER, convert.z.score)
		P_IR<-sapply(Z_IR, convert.z.score)
		P_DR<-sapply(Z_DR, convert.z.score)
	
		mat<-cbind(P_ER, P_IR, P_DR)
		#print(head(mat))
		mat2<-matrix(p.adjust(as.vector(mat), method="bonferroni"),ncol=3)
		adjustedP_ER<-mat2[,1]
		adjustedP_IR<-mat2[,2]
		adjustedP_DR<-mat2[,3]
		
		#ER
		if (sum(is.na(mat2[,1]))) {
		p1ER<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2ER<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("ER",file="Zscore_summary.txt",append=TRUE)
		p1ER<-myplot(Spacings, Z_ER, P_ER, adjustedP_ER, FC_ER)
		p2ER<-volcanoPlot(Spacings, FC_ER, P_ER, adjustedP_ER)
		}
		#IR
		if (sum(is.na(mat2[,2]))) {
		p1IR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2IR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("IR",file="Zscore_summary.txt",append=TRUE)
		p1IR<-myplot(Spacings, Z_IR, P_IR, adjustedP_IR, FC_IR)
		p2IR<-volcanoPlot(Spacings, FC_IR, P_IR, adjustedP_IR)
		}
		#DR
		if (sum(is.na(mat2[,3]))) {
		p1DR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		p2DR<-ggplot() + theme_void() + geom_text(aes(0,0,label='N/A')) + xlab(NULL)
		}else
		{
		write("DR",file="Zscore_summary.txt",append=TRUE)
		p1DR<-myplot(Spacings, Z_DR, P_DR, adjustedP_DR, FC_DR)
		p2DR<-volcanoPlot(Spacings, FC_DR, P_DR, adjustedP_DR)
		}
		pdf(outfig)
# 		print(plot_grid(p1ER, p2ER, p1IR, p2IR, p1DR, p2DR, ncol=2, nrow=3, labels = c('ER', SCO, 'IR', '', 'DR', ''), label_size = 12))
		print(ggdraw()+draw_plot(p1ER, 0.0, 0.667,0.5, 0.333)+draw_plot(p1IR, 0.0, 0.334,0.5, 0.333)+draw_plot(p1DR, 0.0, 0.0,0.5, 0.333)+draw_plot(p2ER, 0.5, 0.667,0.5, 0.333)+draw_plot(p2IR, 0.5, 0.334,0.5, 0.333)+draw_plot(p2DR, 0.5, 0.0,0.5, 0.333)+draw_plot_label(c('ER', 'IR', 'DR', paste("th:",SCO,sep=" ")),c(0.0, 0.0, 0.0, 0.44),c(0.99, 0.67, 0.34, 0.99),size=12))
		dev.off()
# 		+draw_plot_label(c('ER','IR','DR',paste('th:',SCO,sep=" "),c(0.0,0.0,0.0,0.48),c(0.98,0.65,0.0,0.98),size=12))
		fi<-cbind(Spacings, FC_ER,Z_ER,P_ER,adjustedP_ER, FC_IR,Z_IR,P_IR,adjustedP_IR, FC_DR,Z_DR,P_DR,adjustedP_DR)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)

	
	} else if (opt$mat_type=="SYMMETRIC") {
		
		FC_AR<-sapply(Spacings, foldchange, DF=filtDF, CONF="AR", S=Spacings)
		Z_AR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="AR", S=Spacings)
		P_AR<-sapply(Z_AR, convert.z.score)
		adjustedP_AR<-p.adjust(P_AR, method = "bonferroni", n = length(P_AR))
	
		p1<-myplot(Spacings, Z_AR, P_AR, adjustedP_AR, FC_AR)
		p2<-volcanoPlot(Spacings, FC_AR, P_AR, adjustedP_AR)
	
		pdf(outfig, height=4, width=12)
		print(plot_grid(p1, p2, ncol=2, labels = c('A', 'B'), label_size = 12))
		dev.off()
		
		fi<-cbind(Spacings, FC_AR,Z_AR,P_AR,adjustedP_AR)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)
	
	} # AR: All Repeats
	
	outfig2=paste("Zscore_", F, ".png", sep="")
	cmd=paste("convert -density 144", outfig, outfig2, sep=" ")
	system(cmd)
}

print("- - - - - do the job - - - - -")

#infile look like this:
#Peak	Spacing	Score1	Score2	matricePosition1	matricePosition2	correctedPosition1	correctedPosition2
#chr1:19257041-19257441	DR_6	-4.49980967033	-5.19295685089	6	23	17	23
#chr1:19257041-19257441	DR_11	-4.49980967033	-5.29831736655	6	28	17	28
d<-read.table(opt$file, header=TRUE)

#formatage (Jeremy replaced it by a bash command)
# d2 <- d %>% 
# 	separate(Peak, c("tmp", "End"), sep="-") %>%
# 	separate(tmp, c("chr", "Start"), sep=":") %>%
# 	separate(Spacing,c("Conf", "Space"), sep="_")%>%
# 	mutate(Size=as.integer(End)-as.integer(Start)+1)

tx="Listing of over or under represented Spacings for selected scores thresholds"
write(tx,file="Zscore_summary.txt",append=TRUE)

F=1
for (thres in opt$thresholds[[1]]){
	main(d,as.numeric(thres),paste("F",F,sep=""))
	F=F+1
	break
}
warnings()
quit(save="no",status=0)


# extract candidates that are detected in two of three Fscore selections (note uniq -d means at least two)
if (mat_type=="ASYMMETRIC") {
system("cat Zscore_summary.txt | grep ER -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/ER/g' > candidates.txt")
system("cat Zscore_summary.txt | grep IR -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/IR/g' >> candidates.txt")
system("cat Zscore_summary.txt | grep DR -A 1 | grep '[1-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/DR/g' >> candidates.txt")
} else if (mat_type=="SYMMETRIC") {
system("cat Zscore_summary.txt | grep '[1-9]' | tr ' ' '\n' | sort |
uniq -d > candidates.txt")
}

print("finito")