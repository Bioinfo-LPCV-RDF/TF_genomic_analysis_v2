#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
setwd(args[1])
getwd()

# library(ggplot2, lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
# options(bitmapType='cairo')
# library(dplyr, lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
# library(cowplot, lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
# theme_set(theme_cowplot())

library("ggplot2")
options(bitmapType='cairo')
library("dplyr")
library("cowplot")
theme_set(theme_cowplot())


d<-read.table(args[2], header=FALSE)
method=args[3]
name=args[4]

# print("number of sequences:")
n<-nrow(d)
# head(d)
# print(n)


truth<-c(rep(1, n), rep(0, n))


label_prediction<-function(x, score_cutoff) {
	
	x[x>=score_cutoff] <- 1
	x[x<=score_cutoff]<-0
	
	return(x)
}

#set range of scores:
obs_scores<-c(d$V1, d$V2)
# print(obs_scores)
scores<-seq(from = round(min(obs_scores)), to = max(obs_scores), by = 0.01)


compute_measures<-function(SCORE) {

	pred<-c(label_prediction(d$V1, SCORE), label_prediction(d$V2, SCORE))
	
	TP<-length(pred[pred==1 & truth==1]) # or TP<-sum(pred & truth)
	TN<-length(pred[pred==0 & truth==0])
	FP<-length(pred[pred==1 & truth==0])
	FN<-length(pred[pred==0 & truth==1])
	
	all_pred=TP+TN+FP+FN
# 	print(all_pred)
	if (all_pred!=n*2) {print("ERROR"); quit()}
	
	#precision=TP/(TP+FP)
	retrieved <- sum(pred) #i.e. TP+FP
	precision= TP / retrieved

	# recall=TP/(TP+FN)
	recall=TP/sum(truth)   # sum(truth) gives all ground_truths
	
	F1 <- TP / ( TP + 0.5*(FP+FN) )
	F2 <- TP / ( TP + 0.2*FP+ 0.8*FN )
	F3 <- TP / ( TP + 0.1*FP+ 0.9*FN )
	
#  	F1 <- 2 * precision * recall / (precision + recall) #OR TP / (TP + 1/2(FP+FN))
#  	F2 <- 5 * precision * recall / (4 * precision + recall)
	
	
	accuracy<-(TP+TN)/(all_pred)
	specificity<-TN/(TN+FP)
	#return(c(retrieved, precision, recall, F1, F2, specificity))
	return(c(retrieved, precision, recall, F1, F2,F3, specificity))
}


mesu<-t(sapply(scores, compute_measures))
tab<-as.data.frame(cbind(scores, mesu))
colnames(tab)<-c("Score", "retrieved", "precision", "recall", "F1", "F2", "F3", "specificity")
write.table(tab,file=paste("table_fscore", method, ".tsv",sep=""), quote=FALSE, sep="\t")

#print(c(tab$Score[which.max(tab$F1)],tab$Score[which.max(tab$F2)]))
outFile=file(paste("Fscores_", name, "_", method, ".txt",sep=""))
tmp=c(tab$Score[which.max(tab$F1)],tab$Score[which.max(tab$F2)],tab$Score[which.max(tab$F3)])
write(tmp, file = outFile)


# plots
pdf(paste("F_fonc_score", method, ".pdf",sep=""))
#png("F_fonc_score.png", units="px", width=1600, height=1600, res=300)
ggplot(tab, aes(x=Score)) + 
  geom_line(aes(y = F1, linetype="F1")) +
  geom_line(aes(y = F2, linetype="F2")) +
  geom_line(aes(y = precision,  linetype="precision")) +
  geom_line(aes(y = recall,  linetype="recall")) +
  geom_line(aes(y = specificity,  linetype="specificity")) +
  scale_linetype_manual("Measure", values=c("solid", "longdash", "twodash", "dotted", "dashed"))+
  ylab("value") + xlab("score") +
  theme(legend.position=c(0.8, 0.7))
  dev.off()



