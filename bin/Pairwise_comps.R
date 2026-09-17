# Script R : Sample Correlation Analysis
# ======================================

## DESCRIPTION
# This R script generates scatter plots to visualize correlations between samples in a dataset. 
# It uses the psych, ggplot2, and cowplot libraries to create pairwise plots with correlation ellipses and histograms.

## RELATED TO pairwize_comparison

## USAGE
# The script is designed to be called from the command line with the following arguments:
# ```bash
# Rscript script_name.R <working_directory> <data_file> <color>
# ```
## ARGUMENTS
# - `<working_directory>`: The directory where the output plots will be saved.
# - `<data_file>`: The path to the input data file (tab-delimited) containing sample data.

## OUTPUT
# - between_samples_pairwise_scatter_plots.png: Pairwise plot with correlations, ellipses, and histograms.
# - between_samples_logScale.png: Simple pairwise plot with a logarithmic scale.

# LOAD LIBRARIES
.libPaths(c("/home/312.6-Flo_Re/312.6.1-Commun/lib/R-4-4/",.libPaths()))
library(plyr)
library(scales)
library(lazyeval)
library(backports)
library(dplyr)
library(cowplot)
library(ggplot2)
options(bitmapType='cairo') #keep this line
library(psych)
# library(GGally) # Pour les pairwise plots
# CONFIG
theme_set(theme_cowplot())

# LOAD ARGUMENTS
args <- commandArgs(TRUE)
print(args)
setwd(args[1])

# LOAD DATA
d<-read.table(args[2], header=TRUE)

# cols in d
# chr     start   end     ARF10FLrep1     ARF10FLrep2     ARF10FLrep3     From
# chr1    3123    3323    282.122982000554        283.06581131131 339.105517898902        ARF10FLrep1

# from this to the quit('no') => TESTING CODE to replace existing code
# d$From<-factor(d$From)

# # From => sample names
# # if "rep" in sample name, then replicates, then colors in grey; else color in black
# # doing pairwize scatter plots of all samples (cols from 4 to ncol(d)-1)

# data_cols <- colnames(d)[4:(ncol(d) - 1)]

# # Créer une colonne pour la couleur (gris pour "rep", noir sinon)
# d$Color <- ifelse(grepl("rep", d$From), "blue", "red")
# d$alpha <- ifelse(grepl("rep", d$From), 0.3, 0.6)
# # Générer les pairwise scatter plots
# pairwise_plot <- ggpairs(
#   d[, data_cols], # Sélectionner uniquement les colonnes de données
#   aes(color = d$Color), # Colorer les points en fonction de la colonne "Color"
#   upper = list(continuous = wrap("cor", size = 4)), # Corrélation dans les plots supérieurs
#   lower = list(continuous = wrap("points", alpha = d$alpha, size = 1)) # Points dans les plots inférieurs
#   # lower = list(continuous = wrap("points", alpha = 0.3, size = 1)) # Points dans les plots inférieurs
# ) +
#   theme_minimal() +
#   scale_x_log10() + # Échelle logarithmique pour l'axe x
#   scale_y_log10() + # Échelle logarithmique pour l'axe y
#   scale_color_manual(values = c("blue", "red"),name= c("replicate","consensus")) +
#   labs(title = "Pairwise Scatter Plots")

# # Sauvegarder le plot
# ggsave("pairwise_scatter_plots_test.png", pairwise_plot, width = 12, height = 12)

# # Afficher le plot
# # print(pairwise_plot)
# png("between_samples_pairwise_scatter_plots_test.png")
# pairs.panels(d[, data_cols], 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#        log="xy",
#              density = TRUE,  # show density plots
#              ellipses = TRUE, # show correlation ellipses
#              col=d$Color
#              )
# dev.off()

# quit('no')



# check replicats in a pairwise scatter plot
tmp<- d %>% select(-c("chr", "start", "end"))
limMax=max(tmp)
png("between_samples_pairwise_scatter_plots.png")
#pairs(tmp, log="xy")

# PLOTTING
pairs.panels(tmp, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
			 log="xy",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             col=args[3]
             )
dev.off()

png("between_samples_logScale.png")
pairs(tmp, lower.panel=NULL, cex=0.1,col=args[3],xlim=c(0,limMax),ylim=c(0,limMax))
dev.off()