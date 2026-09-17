args <- commandArgs(TRUE)
# library(dplyr, lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
# library(dplyr)
print(args)

infile=args[1]
label=args[2] # this is the mode: inPeaks or inLibs
print(label)
WD=args[3]
print(WD)

setwd(WD)
#load the RC per consensus peaks per replicats
#d<-read.table("allReps_RC.txt", header=F)
d<-read.table(infile, header=F)
# print("coucou")
# head(d)
colnames(d)<-c("chr", "start", "end", args[4:length(args)])
# print("head of read count per peak per sample")
# head(d)
lastCol=ncol(d)


if (label=="inPeaks"){
	#Sum of RC in peaks
	effRC<-colSums(d[,c(4:lastCol)])
	print("total number of reads in peaks per sample")
	print(effRC)
} else if (label=="inLibs") {
	#Final (effective) number of reads in libraries (the ones retained by MACS2)
	RClibs<-read.table("tmpTotalTags.txt", header=FALSE)
# 	RClibs<-read.table("tmpFiltTags.txt", header=FALSE)
	effRC<-RClibs$V1
	print("total number of reads in library per sample")
	print(effRC)
} else {
	print("error: need to tell the program if normalization should be 'inPeaks' or 'inLibs'")
	quit()
}

#normalize peaks RC by tot RC in consensus filtered peaks or lib size
# this give Count Per Million
CPM<-as.data.frame(t(t(d[,c(4:lastCol)]) / (effRC/1000000)))
print("head of RPM per peak per sample")
head(CPM)
# do a RPKM normalization
RPKM<-as.data.frame(CPM / ((d[,3]-d[,2])/1000))
print("head of RPKM per peak per sample")
head(RPKM)
print("summary RPKM")
summary(RPKM)

tmp <- cbind(d[,c(1,2,3)],RPKM)
colnames(tmp)<-c("chr", "start", "end", args[4:length(args)])
filename=paste("peaks_perSample_rpkm", label, ".txt", sep="")
write.table(tmp, file=filename, sep="\t", row.names=F, quote=F)
