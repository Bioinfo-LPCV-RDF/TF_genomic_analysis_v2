

rm(list=ls())
library(DiffBind)

args=commandArgs(trailingOnly=TRUE)
Samplefile=args[1]
outdir=args[2]

samples <- read.csv(Samplefile)
dbObj <- dba(sampleSheet=samples)

# dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.count(dbObj,score=DBA_SCORE_RPKM_FOLD, bUseSummarizeOverlaps=TRUE,summits=100)



png(paste(outdir,"/PCA.png",sep=""))
dba.plotPCA(dbObj)
dev.off()

dbObj <- dba.contrast(dbObj, categories=DBA_TISSUE, minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
print(dbObj)

png(paste(outdir,"/covplot.png",sep=""))
dba.plotMA(dbObj, bXY=TRUE)
dev.off()


res <- dba.report(dbObj, method=DBA_ALL_METHODS, contrast = 1, th=1)
out <- as.data.frame(res)
write.table(out, file=paste(outdir,"/allmethods.tsv",sep=""), sep="\t", quote=F, row.names=F)

res <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(res)
write.table(out, file=paste(outdir,"/Deseq2.tsv",sep=""), sep="\t", quote=F, row.names=F)

res <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(res)
write.table(out, file=paste(outdir,"/edgeR.tsv",sep=""), sep="\t", quote=F, row.names=F)

#  Rscript /home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/bin/Diffbind.R /home/312.6-Flo_Re/312.6.1-Commun/ARF-anr/DAP_052022/results/ARF/ARFanalysis/Comparisons/diffbind_test/ARF5mock_ARF5wIAA_samplesheet.csv /home/312.6-Flo_Re/312.6.1-Commun/ARF-anr/DAP_052022/results/ARF/ARFanalysis/Comparisons/diffbind_test