args <- commandArgs(TRUE)
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
print("THIS IS THE NEW ZACING")
#Rscript Zacing.R SEP3AGspe/SEP3AGspe_spacing_pos.tsv 0.4 yes 0.3 0.4 0.5
#Rscript Zacing.R ARF5_spacing_pos.tsv no -4 -6 -7

# - - - set arguments - - -

infile=args[1]
mat_type=args[2] #should be yes or no

# extacting length of thresholds arrays 
l_thtfa=as.numeric(args[length(args)-1])
l_thtfb=as.numeric(args[length(args)])



# extracting thresholds in lists 
thtfa_list=c()
thtfb_list=c()
for (i in 1:l_thtfa){
	thtfa_list <- append(thtfa_list, as.numeric(args[2+i]))
}
for (i in 1:l_thtfb){
	thtfb_list <- append(thtfb_list, as.numeric(args[l_thtfa+2+i]))
}

FCthreshold=0.6 # 0.6 correspond to FC=1.5
outdir=args[length(args)-2]

# print(args)
setwd(outdir)
meth="median"
if (file.exists("Zscore_summary.txt")){
	system("rm Zscore_summary.txt")
}

# - - load libraries - - -
library(dplyr,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
# library(tidyverse)
library(cowplot,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
theme_set(theme_cowplot())
library(ggrepel,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)

# - - - define functions - - -
foldchange<-function(BP, CONF, DF, S) {

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
	#print("fold change")
# 	print(FC)
	return(FC)	
}

population<-function(BP, CONF, DF, S){
	if (CONF=="AR"){
		DF<-DF
	} else {
		DF<-DF %>% filter(Conf==CONF)
	}
	x<-nrow(DF%>%filter(Space==BP))
	return(x)
}

Zscore_calc<-function(BP, CONF, DF, S) {
#formula adpated from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142402/
# but finaly not really adapted from the article above...;
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
		#print("sigma")
		P=1/length(S)
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


myplot<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,ordmax,ordmin, colourP="grey40", colourL="black") {
	D<-data.frame(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,colourP, colourL)
	
	D$tmp<-NA
	D$tmp[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[abs(log2(D$FC))>FCthreshold & D$ADJUSTED_PVAL<0.05]
	# with the median method I could aslo filter by Zscore> 3.5 as recommended by Iglewicz, B. and Hoaglin, D.C. (2010)
	
	# reporting only positive enrichment
	D$tmp2[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]
	res<-D$tmp2[!is.na(D$tmp2)]
	write(res,file="Zscore_summary.txt",append=TRUE,ncolumns=max(SPA)*2)
	
	p<-ggplot(D, aes(x=SPA, y=ZCO, label=tmp)) +
       geom_point(colour=D$colourP, size=2.5) + geom_line(colour=D$colourL) +
	   geom_text_repel() +
	   ylim(ordmin,ordmax) + 
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}

volcanoPlot<-function(SPA,FC,PVAL, ADJUSTED_PVAL) {
	
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


main<-function(DATA, SCO, SCO2, F) {
# 	print(head(DATA))
	write("- - - - - - - - - - - - - - - - -",file="Zscore_summary.txt",append=TRUE)
	txt=paste(F, SCO, SCO2, sep="_")
	write(txt,file="Zscore_summary.txt",append=TRUE)

	filtDF<-DATA%>%filter(Score1>=SCO & Score2>=SCO2) 
# 	print(head(filtDF))

	minSpa=min(as.numeric(filtDF$Space))
	maxSpa=max(as.numeric(filtDF$Space))
	print(minSpa)
	print(maxSpa)
	Spacings<-seq(minSpa,maxSpa)
	outfig=paste("Zscore_", F, ".pdf", sep="")
	outtable=paste("Zscore_stats_", F, ".tsv", sep="")
	print(mat_type)
	if (mat_type=="ASYMMETRIC") {
		
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

		POP_ER<-sapply(Spacings, population, DF=filtDF, CONF="ER", S=Spacings)
		POP_IR<-sapply(Spacings, population, DF=filtDF, CONF="IR", S=Spacings)
		POP_DR12<-sapply(Spacings, population, DF=filtDF, CONF="DR12", S=Spacings)
		POP_DR21<-sapply(Spacings, population, DF=filtDF, CONF="DR21", S=Spacings)

		P_ER<-sapply(Z_ER, convert.z.score)
		P_IR<-sapply(Z_IR, convert.z.score)
		P_DR12<-sapply(Z_DR12, convert.z.score)
		P_DR21<-sapply(Z_DR21, convert.z.score)
		
		mat<-cbind(P_ER, P_IR, P_DR12, P_DR21)
# 		print(head(mat))
		mat2<-matrix(p.adjust(as.vector(mat), method="bonferroni"),ncol=4)
		adjustedP_ER<-mat2[,1]
		adjustedP_IR<-mat2[,2]
		adjustedP_DR12<-mat2[,3]
		adjustedP_DR21<-mat2[,4]
# 		print(mat2)
		
		ordmax=max(max(Z_ER),max(Z_IR),max(Z_DR12),max(Z_DR21),3,na.rm=TRUE)
		ordmin=min(min(Z_ER),min(Z_IR),min(Z_DR12),min(Z_DR21),0,na.rm=TRUE)
		print(c(ordmin,ordmax))
		
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
		fi<-cbind(Spacings, FC_ER,Z_ER,P_ER,adjustedP_ER,POP_ER, FC_IR,Z_IR,P_IR,adjustedP_IR,POP_IR, FC_DR12,Z_DR12,P_DR12,adjustedP_DR12,POP_DR12,FC_DR21,Z_DR21,P_DR21,adjustedP_DR21,POP_DR21)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)
		

		outfig2=paste("Zscore_", F, ".png", sep="")
		cmd=paste("convert -density 144", outfig, outfig2, sep=" ")
		system(cmd)
		return(list(data.frame("Spacings"=Spacings, "FC_ER"=FC_ER,"Z_ER"=Z_ER,"P_ER"=P_ER,"Padj_ER"=adjustedP_ER, "FC_IR"=FC_IR,"Z_IR"=Z_IR,"P_IR"=P_IR,"Padj_IR"=adjustedP_IR, "FC_DR12"=FC_DR12,"Z_DR12"=Z_DR12,"P_DR12"=P_DR12,"Padj_DR12"=adjustedP_DR12, "FC_DR21"=FC_DR21,"Z_DR21"=Z_DR21,"P_DR21"=P_DR21,"Padj_DR21"=adjustedP_DR21),ordmax,ordmin))

	
	} else if (mat_type=="SYMMETRIC") {
		FC_AR<-sapply(Spacings, foldchange, DF=filtDF, CONF="AR",Spacings)
		Z_AR<-sapply(Spacings, Zscore_calc, DF=filtDF, CONF="AR",Spacings)
		P_AR<-sapply(Z_AR, convert.z.score)
		POP_AR<-sapply(Spacings, population, DF=filtDF, CONF="AR", S=Spacings)
		adjustedP_AR<-p.adjust(P_AR, method = "bonferroni", n = length(P_AR))
	
		ordmax=max(max(Z_AR),3,na.rm=TRUE)
		ordmin=min(min(Z_AR),0,na.rm=TRUE)


		p1<-myplot(Spacings, Z_AR, P_AR, adjustedP_AR, FC_AR,ordmax,ordmin)
		p2<-volcanoPlot(Spacings, FC_AR, P_AR, adjustedP_AR)
	
		pdf(outfig, height=4, width=12)
		print(plot_grid(p1, p2, ncol=2, labels = c('A', 'B'), label_size = 12))
		dev.off()
		
		fi<-cbind(Spacings, FC_AR,Z_AR,P_AR,adjustedP_AR,POP_AR)
		write.table(fi, file=outtable, quote=FALSE, sep="\t", row.names = FALSE)

		outfig2=paste("Zscore_", F, ".png", sep="")
		cmd=paste("convert -density 144", outfig, outfig2, sep=" ")
		system(cmd)
		return(list(data.frame("Spacings"=Spacings, "FC_AR"=FC_AR,"Z_AR"=Z_AR,"P_AR"=P_AR,"Padj_AR"=adjustedP_AR),ordmax,ordmin))
	} # AR: All Repeats
}


myplotg<-function(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC,ordmax,ordmin, Thresholds="black") {
	D<-data.frame(SPA, ZCO, PVAL,ADJUSTED_PVAL, FC, Thresholds)
	
	
	D$ADJUSTED_PVAL <- ifelse(is.finite(D$ZCO), D$ADJUSTED_PVAL,1)
	D$FC <- ifelse(is.finite(D$ZCO), D$FC,1)
	D$ZCO <- ifelse(is.finite(D$ZCO), D$ZCO,0); D$ZCO[is.na(D$ZCO)]<-0
	D$FC <- ifelse(is.finite(D$FC), D$FC,0); D$FC[is.na(D$FC)]<-0
	D$ADJUSTED_PVAL <- ifelse(is.finite(D$ADJUSTED_PVAL), D$ADJUSTED_PVAL,1); D$ADJUSTED_PVAL[is.na(D$ADJUSTED_PVAL)]<-1
	
	D$tmp<-NA
	D$tmp[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]<-D$SPA[log2(D$FC)>FCthreshold & D$ADJUSTED_PVAL<0.05]
	# with the median method I could aslo filter by Zscore> 3.5 as recommended by Iglewicz, B. and Hoaglin, D.C. (2010)
	# print((D))
	p<-ggplot(D,aes(x=SPA,y=ZCO, color=Thresholds,label=tmp)) +
       geom_point(size=2.5) + geom_line() +
	   geom_text_repel() +
	   ylim(ordmin,ordmax) + 
	   xlab("Distance (bp)") + ylab("Z-score")
	return(p)
}

print("- - - - - Zacing for 2 TFBS- - - - -")

#infile look like this:
#Peak	Spacing	Score1	Score2	matricePosition1	matricePosition2	correctedPosition1	correctedPosition2
#chr1:19257041-19257441	DR_6	-4.49980967033	-5.19295685089	6	23	17	23
#chr1:19257041-19257441	DR_11	-4.49980967033	-5.29831736655	6	28	17	28

d<-read.table(infile, header=TRUE)
# print(nrow(d))

# print(head(d))
#formatage
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

for (i in 1:length(thtfa_list)) {
  # create a returnF variable for every threshold and pass them into the main function (to give as 3 dfs)
  returnF_var <- paste("returnF", i, sep="")
  assign(returnF_var, main(d, thtfa_list[[i]], thtfb_list[[i]], paste("F", i, sep="")))

  # create a list_F variable for every threshold to extract the first df in correspondng returnF variable
  list_F_var <- paste("list_F", i, sep="")
  assign(list_F_var, get(returnF_var)[[1]])

  # adding a Thresholds column in df
  list_F_var_df <- get(list_F_var)  # Extract the data frame
  list_F_var_df$Thresholds <- paste("TFa: ",round(thtfa_list[[i]], 2), "; TFb: ",round(thtfb_list[[i]], 2))
  assign(list_F_var, list_F_var_df)

  # chosing ordmax and min (2nd and 3rd dfs in returnF variables)
  ordmax_F_var <- paste("ordmax_F", i, sep="")
  ordmax_F_var <- get(returnF_var)[[2]]
  ordmax <- append(ordmax, ordmax_F_var)
  ordmin_F_var <- paste("ordmin_F", i, sep="")
  ordmin_F_var <- get(returnF_var)[[3]]
  ordmin <- append(ordmin, ordmin_F_var)
	}

# chosing maximum and minimum ord or the different threshold computated dfs
ordmax <- max(ordmax)
ordmin <- min(ordmin)

# merging every thing into a new DF
DF <- do.call(rbind, mget(paste0("list_F", 1:length(thtfa_list))))
DF$Thresholds <- as.factor(DF$Thresholds)

# print(DF)
if (mat_type=="ASYMMETRIC") {
#ER
p1ER<-myplotg(DF$Spacings, DF$Z_ER, DF$P_ER, DF$Padj_ER, DF$FC_ER,ordmax,ordmin, DF$Thresholds)
print("ER")
#IR
p1IR<-myplotg(DF$Spacings, DF$Z_IR, DF$P_IR, DF$Padj_IR, DF$FC_IR,ordmax,ordmin, DF$Thresholds)
print("IR")
#DR
p1DR12<-myplotg(DF$Spacings, DF$Z_DR12, DF$P_DR12, DF$Padj_DR12, DF$FC_DR12,ordmax,ordmin, DF$Thresholds)
print("DR12")
#DR
p1DR21<-myplotg(DF$Spacings, DF$Z_DR21, DF$P_DR21, DF$Padj_DR21, DF$FC_DR21,ordmax,ordmin, DF$Thresholds)
print("DR21")

pdf("Zscore_allF.pdf")
print(plot_grid(p1ER, p1IR, p1DR12, p1DR21, ncol=1, nrow=4, labels = c('ER', 'IR', 'DRab', 'DRba', ''), label_size = 12))
dev.off()

} else if (mat_type=="SYMMETRIC") {

p1AR<-myplotg(DF$Spacings, DF$Z_AR, DF$P_AR, DF$Padj_AR, DF$FC_AR,ordmax,ordmin, DF$Thresholds)
print("AR")

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
if (mat_type=="ASYMMETRIC") {
system("cat Zscore_summary.txt | grep ER -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/ER/g' > candidates.txt")
system("cat Zscore_summary.txt | grep IR -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/IR/g' >> candidates.txt")
system("cat Zscore_summary.txt | grep DRab -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/DRab/g' >> candidates.txt")
system("cat Zscore_summary.txt | grep DRba -A 1 | grep '[0-9]' | tr ' ' '\n' | sort | uniq -d | sed 's/^/DRba/g' >> candidates.txt")
} else if (mat_type=="SYMMETRIC") {
system("cat Zscore_summary.txt | grep '[0-9]' | tr ' ' '\n' | sort |
uniq -d > candidates.txt")
}
# warnings()