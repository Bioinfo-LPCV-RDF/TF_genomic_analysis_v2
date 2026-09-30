.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library("argparse")

parser <- ArgumentParser()

# Rscript /home/312.6-Flo_Re/312.6.1-Commun/ARF-anr/DAP_052022/FeatureAnnotation.R -t1 $outdir/${ARF}mock/${ARF}mock_feature_annotation.tsv -t2 $outdir/${ARF}wIAA/${ARF}wIAA_feature_annotation.tsv -n $ARF -od $outdir/${ARF}/ 

# AIM: plot the distribution of peaks in different genomic features (promoters, exons, introns, intergenic) for the 2 conditions (mock and wIAA) for each ARF

parser$add_argument("-t", "--tableName")
parser$add_argument("-n", "--name")
parser$add_argument("-od", "--outdir")
parser$add_argument("-c", "--color", nargs="+")

args <- parser$parse_args()
tableName <- args$tableName
name <- args$name
out_dir <- args$outdir
color <- args$color

#Creating directory if need be
dir.create(file.path(out_dir),recursive=TRUE ,showWarnings = FALSE)

table <- read.table(tableName,sep="\t",header=FALSE,check.names=FALSE)
colnames(table) <- c("chr","start","end","feature")

table$feature <- factor(table$feature)

# compute the percentage of peaks in each feature for each condition
table <- table %>%
  group_by(feature) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)



# plot the percentages of peaks in different genomic features for the 2 conditions (mock and wIAA) for each ARF
Fontsize=6
p<- ggplot(table, aes(y=feature, x=percentage)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(y="Genomic feature", x="Percentage of peaks") +
  xlim(0, 100) +
  theme_cowplot() +
  theme(
	plot.background = element_rect(fill = "white"),
	axis.text = element_text(size=Fontsize, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "cm")),
	axis.title = element_text(size=Fontsize),
	legend.text = element_text(size=Fontsize),
	legend.title = element_blank(),
	axis.line = element_line(color = "black", linewidth = 0.2),
	axis.ticks = element_line(color = "black", linewidth = 0.2)	
  )  
ggsave(paste0(out_dir, name, "_feature_annotation.png"), p, width=8.4, height=4.2, units="cm",dpi=300)