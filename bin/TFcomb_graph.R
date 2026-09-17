
args <- commandArgs(TRUE)

.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(dplyr)
library(ggplot2)
library(cowplot)
# library(ggnewscale)

# interdistances <- read.table(file="/home/312.6-Flo_Re/312.6.1-Commun/Alice/TFcomb_function/test_from_atlas3/interdistances.txt", sep="\t", header = TRUE)
interdistances <- read.table(file=args[1], sep="\t", header = TRUE)
head(interdistances)

plot <- ggplot(data = interdistances, aes(x = Distance, y = TF1, color=Peak.Heights)) +
        geom_point() +  
        geom_text(aes(label = TF2), hjust = -0.2, vjust = -0.5, size = 3, angle= 70) 
ggsave(paste0(args[2], "TFcomb_interdistances.png"), plot = plot, width = 10, height = 3, bg = 'white')
