rm(list=ls())
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
# library(Biostrings ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library("cowplot",verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(ggplot2 ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
# library(Cairo ,verbose=FALSE ,warn.conflicts=FALSE, quietly = TRUE)
library(RColorBrewer)

args=commandArgs(trailingOnly=TRUE)
TableName=args[1]
out_dir=args[2]
name1=args[3]
name2=args[4]
TableName2=args[5]
suffix=args[6]
maxidiv=args[7]

table<-read.table(TableName2,sep="\t",header=FALSE)
table$V7<-as.factor(table$V7)

Perc<-c()
Name<-c()
Deci<-c()
for (x in 1:maxidiv) { # Percentage according to MACS2
	a<-sum(table$V8[table$V7==x]==name1)
	b<-sum(table$V8[table$V7==x]==name2)
	ab<-sum(table$V8[table$V7==x]=="both")
	aPerc<-(a)/(a+b+ab)
	bPerc<-(b)/(a+b+ab)
	abPerc<-(ab)/(a+b+ab)
	
	Perc<-c(Perc,aPerc*100,abPerc*100,bPerc*100)
	Name<-c(Name,name1,"both",name2)
	Deci<-c(Deci,x,x,x)
}

IntervalsCFC<-c()
Intervals_x<-c()
Intervals_y<-c()
Decil<-c()
meanCFC<-c()
x_vec<-c()
for (x in 1:maxidiv) { 
	miniDec<-round(min(table$V6[table$V7==x]),3)
	maxiDec<-round(max(table$V6[table$V7==x]),3)
	IntervalsCFC<-c(IntervalsCFC,paste("[",miniDec,"-",maxiDec,"]",sep=""))
	Intervals_x<-c(Intervals_x,0.75)
	Intervals_y<-c(Intervals_y,0.024*x+0.05)
	meanDec="CFC ~ 1"
	if(miniDec>1){meanDec="CFC > 1"}
	if(maxiDec<1){meanDec="CFC < 1"}
	
	meanCFC<-c(meanCFC,meanDec)
	Decil<-c(Decil,x)
	x_vec<-c(x_vec,1)
}

MeanDF<-data.frame(as.factor(Decil),as.factor(x_vec),as.factor(meanCFC))
colnames(MeanDF)<-c("Decile","Xval","Attribution")
print(MeanDF)


test_df2<-data.frame(as.factor(Deci),as.factor(Name),as.numeric(Perc))
colnames(test_df2)<-c("Decile","name","Percentage")
colours <- ifelse(as.numeric(test_df2$Percentage) < 50 , "black", "white")

Dataset<-read.table(TableName,sep="\t",header=TRUE)
Dataset$Decile<-as.factor(Dataset$Decile)
Dataset$Spacing<-as.factor(Dataset$Spacing)
print(head(Dataset))
Vmax=max(Dataset$ER,Dataset$DR, Dataset$IR)
Vmin=min(Dataset$ER,Dataset$DR, Dataset$IR)

TrueVmin=max(0,Vmin)
Dataset$ER[Dataset$ER<TrueVmin]=TrueVmin
Dataset$IR[Dataset$IR<TrueVmin]=TrueVmin
Dataset$DR[Dataset$DR<TrueVmin]=TrueVmin

heat1 <- ggplot(Dataset, aes(Spacing, Decile, fill=ER)) + 
geom_tile() + 
ylab(paste(name2," "," ",name1,sep="\t")) + 
scale_fill_gradientn(colors=brewer.pal(n = 9, name = "Reds"),limits=c(TrueVmin,Vmax), name="Z-score")+
theme(panel.background=element_rect(fill = "white"))

heatb<-ggplot(MeanDF, aes(Xval,Decile, fill=Attribution))+
geom_tile()+
scale_color_manual(values = c("CFC ~ 1"="#5bcfcb", "CFC > 1"="#cfb05b", "CFC < 1"="#cf5b97"))+
theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),panel.background=element_rect(fill = "white"),legend.position='none')

heatc<-ggplot(MeanDF, aes(Xval,Decile, fill=Attribution))+
geom_tile()+
scale_color_manual(values = c("CFC ~ 1"="#5bcfcb", "CFC > 1"="#cfb05b", "CFC < 1"="#cf5b97"))+
theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),panel.background=element_rect(fill = "white"))

# scale_fill_gradientn(colors=brewer.pal(n = 3, name = "RdBu"),limits=c(TrueVmin,Vmax), name="Z-score")


# heatb <- ggplot(test_df2, aes(name, Decile, fill=Percentage)) + 
# geom_tile(show.legend=FALSE) + 
# geom_text(aes( label=paste(round(Percentage, 1),"%",sep=""), colour=colours),data=test_df2,size=4,show.legend=FALSE) + 
# xlab("")+
# scale_colour_manual(values = c("white" = "#ffffff", "black" = "#000000"))+
# scale_fill_gradientn(colors=brewer.pal(n = 9, name = "Blues"),limits=c(0,100))+
# theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.background=element_rect(fill = "white"),axis.text.x = element_text(angle = 33, vjust = 1, hjust=1))


heat2 <- ggplot(Dataset, aes(Spacing, Decile, fill=IR)) + 
geom_tile() + 
ylab(paste(name2," "," ",name1,sep="\t")) + 
scale_fill_gradientn(colors=brewer.pal(n = 9, name = "Reds"),limits=c(TrueVmin,Vmax))+
theme(legend.position='none', panel.background=element_rect(fill = "white"))

heat3 <- ggplot(Dataset, aes(Spacing, Decile, fill=DR)) + 
geom_tile() + 
ylab(paste(name2," "," ",name1,sep="\t")) + 
scale_fill_gradientn(colors=brewer.pal(n = 9, name = "Reds"),limits=c(TrueVmin,Vmax))+
theme(legend.position='none', panel.background=element_rect(fill = "white"))

legend1 <-get_legend(heat1)
heat1 <- heat1 + theme(legend.position='none')

legend2 <-get_legend(heatc)


# # GenPlot <- plot_grid(heat1,heatb, heat2,heatb, heat3,heatb, nrow=3)
GenPlot <-ggdraw() +
draw_plot(legend1, 0.85, 0.6, 0.15, 0.5) +
draw_plot(legend2, 0.85, 0.1, 0.15, 0.5) +
draw_plot(heatb, 0.78, 0.035, 0.04, 0.265)+
draw_plot(heatb, 0.78, 0.365, 0.04, 0.265)+
draw_plot(heatb, 0.78, 0.695, 0.04, 0.265)+
draw_plot(heat1, 0.0, 0.66, 0.8, 0.3) +
draw_plot(heat2, 0.0, 0.33, 0.8, 0.3) +
draw_plot(heat3, 0.0, 0, 0.8, 0.3) + 


draw_plot_label(c("DR", "IR", "ER"), c(0.05,0.05,0.05), c(0.31,0.31+0.33,0.31+0.66), size = 12)

# draw_plot_label(IntervalsCFC, Intervals_x, Intervals_y, size = 8)+
# draw_plot_label(IntervalsCFC, Intervals_x, Intervals_y+0.33, size = 8)+
# draw_plot_label(IntervalsCFC, Intervals_x, Intervals_y+0.66, size = 8)+
# draw_plot_label(c("CFC intervals","CFC intervals","CFC intervals"), c(0.74,0.74,0.74), c(0.31,0.31+0.33,0.31+0.66), size = 9)+




save_plot(paste(out_dir,"/Heatmap_decile_",suffix,".png",sep=""), GenPlot, base_width=8, base_height=9
          )

GenPlot3 <- plot_grid(heat1, heat2,heat3, nrow=3)

save_plot(paste(out_dir,"/Heatmap_decile_",suffix,"oldformat.png",sep=""), GenPlot3, base_width=8, base_height=9
          )
          
table<-read.table(TableName2,sep="\t",header=FALSE)
table$V7<-as.factor(table$V7)

cov1<-ggplot(data=table,aes(y=V4,x=V5))+
	geom_point(alpha=0.4,aes(color=V7)) +
	scale_x_log10() +
	scale_y_log10() +
	labs(color=paste("Decile of peaks\n with ratio\n ",name1,"/",name2,sep="")) +
	# scale_color_manual(values=c("#8243D6","#39E639","#FF9640"))+
	# scale_color_manual(values = c(paste(name1) = "#7fc3c0", paste(name2) = "#cfb845", "both" = "#141414"))+
# 	theme(axis.title=element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=12)) +
	# guides(colour = guide_legend(override.aes = list(size=10,alpha=1))) +
	geom_abline(intercept = 0) +
	xlab (paste("cov ",name2,"",sep="")) +
	ylab (paste("cov ",name1,"",sep=""))
GenPlot2<- plot_grid(cov1)

save_plot(paste(out_dir,"/cov_decile_",suffix,".png",sep=""), GenPlot2, base_width=8, base_height=6
          )