#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library(ggplot2)
library(Cairo)
library(dplyr)
library(cowplot, lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
theme_set(theme_cowplot())

ratio1=as.numeric(args[1])
print(ratio1)
ratio2=as.numeric(args[2])
print(ratio2)
label=args[3]
print(label)
WD=args[4]
print(WD)
numRep1=as.integer(args[5])
print(numRep1)
numRep2=as.integer(args[6])
print(numRep1)
nameCond1=args[7]
nameCond2=args[8]
print("replicats names:")
print(args[9:14])

# in this script we ll say that:
# M for mutant is SEP3AGiAP1
# W for wilde is SEP3AG

setwd(WD)
#load the RC per consensus peaks per replicats
d<-read.table("allReps_RC.txt", header=F)
head(d)
head(d[,2])
#colnames(d)<-c("chr", "start", "end", "W1", "W2", "W3", "M1", "M2", "M3")
#colnames(d)<-c("chr", "start", "end", "D1", "D2", "A1", "A2", "A3")
#summary(d$end-d$start)
#summary(d[,c(4:9)])

if (label=="inPeaks"){
#Sum of RC in peaks
lastCol=4+numRep1+numRep2-1
effRC<-colSums(d[,c(4:lastCol)])
print(effRC)
} else if (label=="inLibs") {
	#Final (effective) number of reads in libraries (the ones retained by MACS2)
	RClibs<-read.table("tmpFiltTags.txt", header=FALSE)
	effRC<-RClibs$V1
} else {
	print("error: need to tell the program if normalization should be 'inPeaks' or 'inLibs'")
	quit()
}

#normalize peaks RC by tot RC in consensus filtered peaks or lib size
# this give Count Per Million
CPM<-as.data.frame(t(t(d[,c(4:lastCol)]) / (effRC/1000000)))
head(CPM)
# do a RPKM normalization
RPKM<-as.data.frame(CPM / ((d[,3]-d[,2])/1000))
#print("summary RPKM")
head(RPKM)
head(RPKM[,1])

tmp <- cbind(d[,c(1,2,3)],RPKM)
colnames(tmp)<-c("chr", "start", "end", args[9:14])
write.table(tmp, "mergedPeaks_perSampRPKMRIP.txt", sep="\t", row.names=F, quote=F)


covPlot<-function(rep1, rep2, YNAME, XNAME) {
	d1<-data.frame(rep1, rep2)
	colnames(d1)<-c("rep1", "rep2")
	d2<- data.frame(d1 %>% mutate(catego=cut(rep1/rep2, 
		breaks=c(-Inf, ratio1, ratio2 , Inf),
		labels=c("rep2","Shared","rep1"))))
	#print(d2%>%filter(catego!="Shared"))
	TS=10
	
	t=cor.test(d2$rep1, d2$rep2, method="pearson", use = "complete.obs")
	r=round(t$estimate, 2)
	P=round(t$p.value, 10)
	#P=round(log10(t$p.value), 2)
	
	p <- ggplot(data=d2 ,aes(y=rep1,x=rep2, colour=factor(catego))) +
		scale_x_log10() + 
		scale_y_log10() + 
		xlab(bquote(paste(.(XNAME), ' coverage ( log'['10']*'(RPKM) )'))) + 
		ylab(bquote(paste(.(YNAME), ' coverage ( log'['10']*'(RPKM) )'))) +
		#xlab(paste(XNAME, 'coverage (RPKM)')) +
		#ylab(paste(YNAME, 'coverage (RPKM)')) +
		theme(legend.position="none",axis.text=element_text(size=TS),
				axis.title=element_text(size=TS)) +
		guides(color = guide_legend(override.aes = list(size = 4))) +		
		scale_color_manual(name="", values = c("rep1" = "#d7191c",
			"Shared" = "#d9d9d9", "rep2"="#d7191c")) +
		geom_point(size=0.3) +
		geom_abline(intercept = 0, linetype = "dashed", color="darkgrey") +
		#scale_x_log10(breaks=c(1,10,100),labels=c(1, 10, 100)) +
		#scale_y_log10(breaks=c(1,10,100),labels=c(1, 10, 100)) +
		##2500 et 3000 en x et y lim pour in Libs
		#xlim(0,6000) +
		#ylim(0,6000) +
		#annotate('text', 1000, 2500,
		#	label=paste("italic(r)==", r, "~italic(log10(P))==", P), parse=TRUE,
		#	hjust=1, size=4)
		annotate('text', 1000, 5500,
			label=paste("italic(r)==", r), parse=TRUE,
			hjust=1, size=4)
	
	return(p)
}


if (numRep1==3 && numRep2==3) {
	print("hello")
	p1<-covPlot(RPKM[,1], RPKM[,2], args[9], args[10])
	p2<-covPlot(RPKM[,1], RPKM[,3], args[9], args[11])
	p3<-covPlot(RPKM[,2], RPKM[,3], args[10], args[11])
	
	p4<-covPlot(RPKM[,4], RPKM[,5], args[12], args[13])
	p5<-covPlot(RPKM[,4], RPKM[,6], args[12], args[14])
	p6<-covPlot(RPKM[,5], RPKM[,6], args[13], args[14])
	
	pdf(paste("pairwise_rep_coverage_", label, ".pdf", sep=""), width=16, height=8)
    plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D', 'E', 'F'))
	#dev.off()
} else if (numRep1==2 & numRep2==3) {
	p1<-covPlot(RPKM[,4], RPKM[,5], args[9], args[10])
	
	p1<-covPlot(RPKM[,6], RPKM[,7], args[11], args[12])
	p2<-covPlot(RPKM[,6], RPKM[,8], args[11], args[13])
	p3<-covPlot(RPKM[,7], RPKM[,8], args[12], args[13])
	
	pdf(paste("pairwise_rep_coverage_", label, ".pdf", sep=""), width=10, height=8)
    plot_grid(p1, p4, p5, p6,labels = c('A', 'B', 'C', 'D'))
    dev.off()
} else {print("this code is needs to be arrange for your number of replicates per condition")}
#system("convert -density 144 pairwise_rep_coverage_inLibs.pdf pairwise_rep_coverage_inLibs.png")
#system("convert -density 144 pairwise_rep_coverage_inPeaks.pdf pairwise_rep_coverage_inPeaks.png")


# conditions plots
#tmp<-data.frame(rowMeans(RPKM[c("V4","V5","V6")]), 
	#rowMeans(RPKM[c("V7","V8","V9")]))
tmp<-data.frame(rowMeans(RPKM[,c(1:numRep1)]), 
	rowMeans(RPKM[,c((numRep1+1):(numRep1+numRep2))]))
colnames(tmp)<-c("exp1", "exp2")
head(tmp)

d2<- data.frame(tmp %>% mutate(catego=cut(exp1/exp2, 
		breaks=c(-Inf, ratio2, ratio1, Inf),
		labels=c("exp2","Shared","exp1"))))
print(nrow(d2%>%filter(catego=="exp2")))
head(d2)
summary(d2)
nrow(d2)

TS=14

p <- ggplot(data=d2 ,aes(y=exp1,x=exp2, colour=factor(catego))) +
	scale_x_log10() + 
	scale_y_log10() + 
	xlab(bquote(paste(.(nameCond2), ' coverage ( log'['10']*'(RPKM) )'))) + 
	ylab(bquote(paste(.(nameCond1), ' coverage ( log'['10']*'(RPKM) )'))) +
	#xlab(expression("coverage ( "*log[10]~"(RPKM) )")) + 
	#ylab(expression("coverage ( "*log[10]~"(RPKM) )")) +
	#theme_bw() +
	theme(legend.position="none",axis.text=element_text(size=TS),
			axis.title=element_text(size=TS)) +
	guides(color = guide_legend(override.aes = list(size = 4))) +		
	scale_color_manual(name="", values = c("exp1" = "#ef8a62",
		"Shared" = "#d9d9d9", "exp2"="#67a9cf")) +
	geom_point(size=0.3) +
	geom_abline(intercept = 0, linetype = "dashed", color="darkgrey")
	#scale_x_log10(breaks=c(1,10,100),labels=c(1, 10, 100))

h<-ggplot(data=d2 ,aes(exp1/exp2, colour=factor(catego), fill=factor(d2$catego))) +
	geom_histogram(binwidth=0.05, alpha=.6, position="identity") +
	scale_x_log10() + 
	#theme_bw() + 
	theme(legend.position=c(.7,.9),axis.text=element_text(size=TS),
			axis.title=element_text(size=TS), legend.background=element_blank(),
			legend.key.width = unit(.5, "cm")) +
	scale_color_manual(name="", values = c(exp1 = "black", "Shared" = "black", exp2="black"), labels=c(nameCond2, "Shared", nameCond1)) +
	scale_fill_manual(name="", values = c(exp1 = "#ef8a62", "Shared" = "#d9d9d9", exp2="#67a9cf"), labels=c(nameCond2, "Shared", nameCond1))
	#scale_fill_discrete(name = "", labels = c("A", "B", "C"))
#	xlab("CFR nameCond1/nameCond2") +
#	ylab("Count")

pdf(paste("conditions_coverage_", label,".pdf", sep=""), width=10, height=5)
plot_grid(p, h, labels = c('A', 'B'))
dev.off()
#system("convert -density 144 conditions_coverage_inLibs.pdf conditions_coverage_inLibs.png")
#system("convert -density 144 conditions_coverage_inPeaks.pdf conditions_coverage_inPeaks.png")

out <- cbind(d[,c(1,2,3)],d2, d2$exp1/d2$exp2)
colnames(out)<-c("chr","start", "end", nameCond1, nameCond2, "category", "CFC")
#out <- out[order(,"CFC", decreasing=FALSE)]
print(head(out))
write.table(out, "avRPKM_RCinPeaks.txt", sep="\t", row.names=F, quote=F)

quit()

# PLOTS WITH ERRORS BARS # needs to be updated
d2<-data.frame(rowMeans(RPKM[c("D1", "D2")]),
        rowMeans(RPKM[c("A1", "A2", "A3")]))
colnames(d2)<-c("SEP3delAG", "SEP3del3mAG")

sdx <- data.frame(apply(as.matrix(RPKM %>% select(c("D1","D2"))),1,FUN= sd))
y1<-d2[,"SEP3delAG"]
sdy <- data.frame(apply(as.matrix(RPKM %>% select(c("A1","A2","A3"))),1,FUN= sd))
x1<-d2[,"SEP3del3mAG"]

colnames(sdx) <- 'sdx'
colnames(sdy) <- 'sdy'
table <- cbind(d[,c(1,2,3)],d2,sdx,sdy)
table <- table[order(sdx/x1+sdy/y1,decreasing=FALSE),]
write.table(table, "all_avRPKM_RCinPeaks.txt", sep="\t", row.names=F, quote=F)

print("summary table:")
summary(table)
commonFold<-table %>% filter(SEP3del3mAG/SEP3delAG<args[2] & SEP3del3mAG/SEP3delAG>args[1]) %>% select(chr, start, end)
write.table(commonFold, sep="\t", file=paste("common",label,".bed",sep=""), row.names=F, quote=F, col.names=F)

table<- data.frame(table %>% mutate(category=cut(SEP3del3mAG/SEP3delAG, breaks=c(-Inf, args[2], args[1], Inf), labels=c("SEP3delAG","Shared","SEP3del3mAG"))))

print("coucou")

pdf(paste("biplot_avSEP3delAG_vs_avSEP3del3mAG_RPKM_errorBar_", label,".pdf", sep=""))
#Cairo(width = 1500, height = 1000, file="biplot_avAG_vs_avAGi_RPKM_errorBar.png", type="png", pointsize=12, bg = "transparent", canvas = "white", units = "px", dpi = "auto")
g <- ggplot(data=table ,aes(x=SEP3delAG,y=SEP3del3mAG, colour=factor(category)))
g + scale_x_log10() + 
	scale_y_log10() + 
	ylab("SEP3del3mAG relative binding intensity") + 
	xlab("SEP3delAG relative binding intensity") +
	#theme(axis.title=element_text(size=20),legend.text=element_text(size=25)) + 
	#guides(colour = guide_legend(override.aes = list(size=10,alpha=1))) + 
    guides(shape=FALSE) + 
	theme(legend.position="none",axis.text=element_text(size=24), axis.title=element_text(size=24)) +
	geom_errorbarh(height=0,aes(xmin=SEP3delAG-sdx, xmax=SEP3delAG+sdx)) +
	geom_errorbar(width=0,aes(ymin=SEP3del3mAG-sdy, ymax=SEP3del3mAG+sdy)) +
	scale_color_manual(name="CFR category", values = c("SEP3del3mAG" = "#000000", "Shared" = "#d9d9d9", "SEP3delAG"="#000000")) +
	geom_point() +
	geom_abline(intercept = 0)
dev.off()
system("convert -density 144 biplot_avSEP3delAG_vs_avSEP3del3mAG_RPKM_errorBar_inPeaks.pdf biplot_avAG_vs_avAGi_RPKM_errorBar_inPeaks.png")
system("convert -density 144 biplot_avSEP3delAG_vs_avSEP3del3mAG_RPKM_errorBar_inLibs.pdf biplot_avSEP3delAG_vs_avSEP3del3mAG_RPKM_errorBar_inLibs.png")

quit()

# - - - - - - - output outliers as a bed - - - - - - -

AG_foldOutliers<- table %>% filter(AGi/AG<args[2]) %>% select(chr, start, end)
write.table(AG_foldOutliers, sep="\t", file=paste("AG_foldOutliers_",label,".bed",sep=""), row.names=F, quote=F, col.names=F)

AGi_foldOutliers<- table %>% filter(AGi/AG>args[1]) %>% select(chr, start, end)
write.table(AGi_foldOutliers, sep="\t", file=paste("AGi_foldOutliers_",label,".bed",sep=""), row.names=F, quote=F, col.names=F)


#sorted
sortedAGitoAG<-table %>% mutate(AGitoAG=AGi/AG) %>% arrange(desc(AGitoAG)) %>% mutate(sdx=round(sdx, 1), sdy=round(sdy, 1))
write.table(sortedAGitoAG, sep="\t", file=paste("sortedAGitoAG_",label,".txt",sep=""), row.names=F, quote=F, col.names=T)

f=0.2

topAG=tail(sortedAGitoAG, nrow(sortedAGitoAG)*f)
write.table(topAG, sep="\t", file=paste("topAG_",label,".txt",sep=""), row.names=F, quote=F, col.names=T)

topAGi=head(sortedAGitoAG, nrow(sortedAGitoAG)*f+1)
write.table(topAGi, sep="\t", file=paste("topAGi_",label,".txt",sep=""), row.names=F, quote=F, col.names=T)

































