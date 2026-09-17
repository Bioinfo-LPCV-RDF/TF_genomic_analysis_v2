print("___  ___           __   __         __  ")
print(" |  |__     __    /  ` /  \  |\/| |__) ")
print(" |  |             \__, \__/  |  | |__) ")
print("                                       ")

# TO DO :
#_______
# Review selection of significant rules 
# add an argument with the name of the TF for the interdistance analysis being flexible
# sortir plot zscore cosine significant rule avec p value (voir docu)
# Interdistance: Why does a peak surpass the Zscore threshold but it is not considered as a significant interdistance? 

"""
TF-COMB is a program predicting co-occurences and interdistances between TF hetero-dimer under peaks

Analysis is splited into 2 main parts: 
    - Co-occurrence calculation -> result hetero-dimer TF pairs named "rules" 
    - Interdistance calculation from previous analysis  

Arguments: 
    - peaks: peaks (narrow peaks from ChIP sample in .bed)
    - motifs_file: concatenation of matrices into a single file (.txt)
    - genome: reference genome (.fa)
    - m_pvalue: set a pvalue to filter matrices 
    - pvalue: set the threshold taken to select significant rules 
    - TH: thread number to perform analysis 

Original GitHub: https://github.com/loosolab/TF-COMB/
Documentation: https://tf-comb.readthedocs.io/en/latest/

Author: Alice JEGOU  
Date: December 2024
"""

# INITIALIZATION: 
# _______________

# Importations
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_rows', None)
from tfcomb import CombObj
C = CombObj()

# Load arguments
peaks=sys.argv[1]
print("peaks file: ", peaks)
motifs_file=sys.argv[2]
print("motifs file: ", motifs_file)
genome=sys.argv[3]
print("genome file: ", genome)
pvalue = float(sys.argv[4])
print("pvalue: ", pvalue)   
m_pvalue = float(sys.argv[5])
th=int(sys.argv[6])
out_dir=Path(sys.argv[7])
TF_name=sys.argv[8]

# Create a TF-COMB object
print("--- Creating TFComb Object ---")
C.TFBS_from_motifs(regions=peaks,
                   motifs=motifs_file,
                   genome=genome,
                   threads=th, motif_pvalue=m_pvalue, 
                   resolve_overlapping='highest_score')


# CO-OCCURRENCE ANALYSIS:
#________________________

print("--- Co-occurence analysis ---")

# List of TFBS found under peaks (strand is precised by default)
output_file = out_dir / 'TFBS.txt'
with open(output_file, 'w') as file:
    file.write(str(C.TFBS))

# Co-occurence calculation
C.count_within(max_overlap=1, threads=th, binarize=True)
C.market_basket(threads=th)
C.simplify_rules() #deduplicate rules (easier to read)

# write output file listing hetero-dimer pairs 
output_file = out_dir / 'co_occurrence_rules.txt'
C.rules.to_csv(output_file, sep='\t', index=False)
# Select significant rules
selection = C.select_significant_rules(x="zscore",y="zscore",x_threshold_percent=pvalue, y_threshold_percent=pvalue)
output_file = out_dir / 'co_occurrence_rules_significant.txt'
selection.rules.to_csv(output_file, sep='\t', index=False)



# _=C.plot_bubble(selection)
# heatmap_figure = plt.gcf()
# if heatmap_figure:
#     heatmap_figure.savefig("figure_bubble.png", format="png")


# INTER-DISTANCES ANALYSIS:
#__________________________

print("--- Interdistances analysis ---")

selection.analyze_distances(threads=th)
# Write out results
output_file = out_dir / 'interdistances.txt'
with open(output_file, 'w') as file:
    selection.distObj.peaks.loc[(selection.distObj.peaks.TF1 == TF_name)].to_csv(output_file, sep='\t', index=False)
    # file.write(str(selection.distObj.peaks.loc[(selection.distObj.peaks.TF1 == TF_name)]))

# _ = selection.distObj.plot(("LFY", "SEP3"))
# output_file = out_dir / 'interdistances.png'
# spacing_figure = plt.gcf()
# if spacing_figure:
#     spacing_figure.savefig(output_file, format="png")



