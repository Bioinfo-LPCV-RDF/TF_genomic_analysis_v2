# coding: utf-8
"""
get_all_score_tffm.py 
======================

This script processes a TFFM (Transcription Factor Flexible Model) file and scans a given 
FASTA file for sequences that match the TFFM model. The results are written to an output file.

Related to: compute_distribution, compute_space, compute_spacing2TFs, compute_spacing_TFBSTSS

Dependencies:
- tffm_module: A module for handling TFFM models.
- constants: Contains the TFFM_KIND constant.

Arguments:
- --output, -o: Path to the output file where results will be written.
- --fasta_pos, -pos: Path to the input FASTA file containing sequences to scan.
- --tffm, -t: Path to the TFFM XML file.
- --threshold, -th: Threshold score for sequence matching (default: 0).

Functionality:
1. Loads a TFFM model from the specified XML file.
2. Scans the sequences in the provided FASTA file using the TFFM model.
3. Writes the positions of matching sequences to the output file.
Usage:
Run the script with the required arguments to perform the sequence scan and save the results.

"""
# IMPORTATIONS
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

parser = argparse.ArgumentParser()                                               

parser.add_argument("--output", "-o", type=str)
parser.add_argument("--fasta_pos", "-pos", type=str)
parser.add_argument("--tffm", "-t", type=str)
parser.add_argument("--threshold", "-th", type=float,default=0)
args = parser.parse_args()
output = args.output
tffm=args.tffm
fasta_pos=args.fasta_pos
threshold=args.threshold

hit_pos_list=[]
tffm_first_order = tffm_module.tffm_from_xml(tffm, TFFM_KIND.FIRST_ORDER)
print (tffm)
print (fasta_pos)

with open(output,"w") as f1:
    for hit_pos in tffm_first_order.scan_sequences(fasta_pos,threshold=threshold, only_best=False):
#        hit_pos_list.append(str(hit_pos))
         if str(hit_pos)!="None":

             f1.write(str(hit_pos)+"\n")        
    
#with open(output,"w") as f1:
#    for hit_pos in hit_pos_list:
#        f1.write(hit_pos+"\n")        
