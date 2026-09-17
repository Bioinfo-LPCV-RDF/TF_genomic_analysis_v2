#!/usr/bin/env python
# coding: utf-8

"""
TFFM Model Training and Visualization
=====================================

DESCRIPTION:
This script trains a TFFM (Transcription Factor Flexible Model) using a MEME motif and FASTA sequences.
It generates both initial and trained TFFM models, and creates visual representations (logos) of these models.

RELATED TO compute_motif

USAGE:
```bash
python tffm_train.py --results <results_dir> --fasta <fasta_file> --matrix <meme_matrix> --place <motif_place>
```

ARGUMENTS:
Argument,Description
- --results, -r, Directory where results will be saved. Required.
- --fasta, -f, Path to the FASTA file containing sequences for training. Required.
- --matrix, -m, Path to the MEME motif file. Required.
- --place, -p, Index of the motif to use from the MEME file. Required.

OUTPUTS: 
- tffm_summary_initial.svg: Summary logo of the initial TFFM model.
- tffm_dense_initial.svg: Dense logo of the initial TFFM model.
- tffm_first_order_initial.xml: Initial TFFM model in XML format.
- tffm_first_order.xml: Trained TFFM model in XML format.
- tffm_summary.svg: Summary logo of the trained TFFM model.
- tffm_dense.svg: Dense logo of the trained TFFM model.

NOTES:
- The MEME motif file should be in standard MEME format.
- The FASTA file should contain sequences relevant to the motif being analyzed.
- The place argument specifies which motif to use if the MEME file contains multiple motifs.
"""

# CONFIG 
from os import system
from os import mkdir
from os.path import isdir
from os.path import basename
from os.path import dirname
from os.path import join
import sys
sys.path.append("/home/312.3-StrucDev/312.3.1-Commun/lib/TFFM")
import tffm_module
from constants import TFFM_KIND
import argparse

# LOAD ARGUMENTS
parser = argparse.ArgumentParser()                                               

parser.add_argument("--results", "-r", type=str)
parser.add_argument("--fasta", "-f", type=str)
parser.add_argument("--matrix", "-m", type=str)
parser.add_argument("--place", "-p", type=int)
args = parser.parse_args()
results = args.results
meme_motif=args.matrix
fasta=args.fasta
place=args.place

# print(fasta)
# print(results)
# print(meme_motif)
# print(place)

# TFFM MODELLING
tffm_first_order = tffm_module.tffm_from_meme(meme_motif, TFFM_KIND.FIRST_ORDER,place)
out = open(results+"/tffm_summary_initial.svg", "w")
tffm_first_order.print_summary_logo(out)
out.close()
tffm_first_order.write(results+"/tffm_first_order_initial.xml")
out = open(results+"/tffm_dense_initial.svg", "w")
tffm_first_order.print_dense_logo(out)
out.close()
tffm_first_order.train(fasta)
tffm_first_order.write(results+"/tffm_first_order.xml")
out = open(results+"/tffm_summary.svg", "w")
tffm_first_order.print_summary_logo(out)
out.close()
out = open(results+"/tffm_dense.svg", "w")
tffm_first_order.print_dense_logo(out)
out.close()









