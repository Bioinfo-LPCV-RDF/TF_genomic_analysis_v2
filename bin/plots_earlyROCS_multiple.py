#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ROC Curve Generator and AUC Calculator
========================================

DESCRIPTION:
This script generates ROC (Receiver Operating Characteristic) curves and calculates
AUC (Area Under the Curve) values for multiple score datasets. It's designed to compare
the performance of different models or methods in distinguishing between bound and
unbound sequences.

The script reads score files containing positive and negative scores, calculates ROC
curves for each dataset, and generates a combined plot with all curves and their AUC values.

AUTHOR: Jérémy LUCAS (inspired by Arnaud STIGLIANI)
DATE: 17-01-2018

USAGE:
The script is designed to be called from the command line with the following arguments:
```bash
python roc_curve_generator.py --scores <score_files> --names <dataset_names>
                              --colors <colors> --output <output_dir>
                              --outputfile <output_filename>
```

ARGUMENTS:
--scores, -s       : List of score files to process. Required.
--names, -n        : List of names for each dataset. Required.
--colors, -c       : List of colors for plotting each curve. Required.
--output, -o       : Output directory path. Required.
--outputfile, -of  : Output filename for the plot. Required.

FUNCTIONS: 
- prepare_curve()  : Calculates FPR, TPR, and AUC for a dataset
- plot_curve()     : Generates the ROC curve plot with all datasets

INPUT FILES: 
Each score file should contain two columns:
- Column 1: Positive sequence scores
- Column 2: Negative sequence scores

Example Input:
0.95    0.12
0.87    0.23
0.92    0.18
...

OUTPUTS: 
- <output_dir>/<output_filename>  : PNG image of the ROC curves
- <output_dir>/AUC.txt            : Text file with AUC values for each dataset
- <output_dir>/nKSM_AUC.txt       : Text file with numerical AUC values

NOTES: 
- The script uses matplotlib's Agg backend for non-interactive plotting.
- Each dataset should have the same number of positive and negative scores.
- The diagonal line (random classifier) is included in the plot for reference.
- AUC values are rounded to 3 decimal places for display.
- The script handles multiple datasets for comparison on the same plot.
"""

################################################                import packages                ################################################                

import argparse
import os
import matplotlib #↑ 2.2.5
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sklearn.metrics # 0.20.4
from sklearn.metrics import roc_curve, auc
import threadpoolctl 
import numpy as np # 1.16.6
################################################                Parsing arguments                ################################################

parser = argparse.ArgumentParser() 

parser.add_argument("--scores","-s" , nargs="+")
parser.add_argument("--names","-n" , nargs="+")
parser.add_argument("--colors","-c" , nargs="+")
parser.add_argument("--output","-o")
parser.add_argument("--outputfile","-of")
parser.add_argument("--fpr_limit","-lim")
args = parser.parse_args()
name_scores_files = args.scores
names = args.names
dir_OUT = args.output
fileout = args.outputfile
colors = args.colors
fpr_limit = float(args.fpr_limit)

################################################                Functions                ################################################'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black']
def prepare_curve(groups, scores, fpr_limit):
    fpr, tpr, _ = roc_curve(groups, scores, pos_label=1)

    # Keep only the early part of the ROC curve (FPR <= fpr_limit)
    mask = fpr <= fpr_limit

    # Interpolate TPR at exact fpr_limit
    fpr_partial = np.concatenate([fpr[mask], [fpr_limit]])
    tpr_partial = np.concatenate([tpr[mask], [np.interp(fpr_limit, fpr, tpr)]])

    early_auc = auc(fpr_partial, tpr_partial)
    early_auc = round(early_auc, 3)

    return fpr, tpr, early_auc

def plot_early_roc(scores, name_scores_files, names, dir_OUT, fileout, colors, fpr_limit):
    fig = plt.figure(figsize=(6, 6), dpi=60)
    ax = fig.add_subplot(111)

    ax.set_xlabel("fraction of unbound sequences", fontsize=14)
    ax.set_ylabel("fraction of bound sequences", fontsize=14)

    fo = open(os.path.join(dir_OUT, "Early_AUC.txt"), "w")

    for i, elt in enumerate(name_scores_files):

        # Filter only on FPR (same thresholds, no TPR filtering)
        mask = scores[elt]["fpr"] <= fpr_limit
        fpr_early = scores[elt]["fpr"][mask]
        tpr_early = scores[elt]["tpr"][mask]

        ax.plot(
            fpr_early,
            tpr_early,
            linewidth=3,
            color=colors[i],
            label=names[i]
        )

        ax.annotate(
            "AUC " + names[i] + " = " + str(scores[elt]["auc"]),
            xy=(0.02, 0.05 + 0.05 * i),
            color=colors[i],
            fontsize=14
        )

        print("AUC " + names[i] + " = " + str(scores[elt]["auc"]))
        fo.write("AUC " + names[i] + " = " + str(scores[elt]["auc"]) + "\n")

    # Random classifier diagonal (early region)
    ax.plot([0, fpr_limit], [0, fpr_limit], linestyle="dashed", color="black")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(dir_OUT, fileout))
    fo.close()
################################################                MAIN                ################################################

scores = {}
fout = open(os.path.join(dir_OUT, "nKSM_AUC.txt"), "w")

for elt in name_scores_files:
    scores[elt] = {
        "groups": [],
        "scores": [],
        "fpr": [],
        "tpr": [],
        "auc": 0
    }

    with open(elt, "r") as IN:
        for line in IN:
            columns = line.strip().split("\t")
            if len(columns) >= 2 and columns[0] not in ["NA", ""] and columns[1] not in ["NA", ""]:
                scores[elt]["groups"].append(1)  # positive
                scores[elt]["scores"].append(float(columns[0]))
                scores[elt]["groups"].append(0)  # negative
                scores[elt]["scores"].append(float(columns[1]))

    fpr, tpr, early_auc = prepare_curve(scores[elt]["groups"], scores[elt]["scores"], fpr_limit)
    

    scores[elt]["fpr"] = fpr
    scores[elt]["tpr"] = tpr
    scores[elt]["auc"] = early_auc

    fout.write("%f\n" % early_auc)

fout.close()

plot_early_roc(scores, name_scores_files, names, dir_OUT, fileout, colors, fpr_limit)

#scores_1=[]
#groups_1=[]
#fpr_1=""
#tpr_1=""
#auc_1=""

#IN.close()

#scores_2=[]
#groups_2=[]
#fpr_2=""
#tpr_2=""
#auc_2=""
#with open(name_scores_file_2,"r") as IN:
#    for line in IN:
#        columns=line.split("\t")
#        if columns[0] != "NA" and columns[1] != "NA":
#            groups_2.append(1)
#            scores_2.append(float(columns[0]))
#            groups_2.append(2)
#            scores_2.append(float(columns[1]))
#IN.close()

#(fpr_1,tpr_1,auc_1)=prepare_curve(groups_1,scores_1)
#(fpr_2,tpr_2,auc_2)=prepare_curve(groups_2,scores_2)
#plot_curve(fpr_1,tpr_1,auc_1,fpr_2,tpr_2,auc_2,name_1,name_2)

