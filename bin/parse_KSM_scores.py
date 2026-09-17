#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
KSM Scores Parser and Comparator
===============================

DESCRIPTION
This script parses KSM (k-mer set motif) score files and compares scores between positive
and negative sequence sets. It extracts scores for each sequence, calculates summary
statistics (sum, max, mean) for each sequence's scores, and outputs a comparison file
with matched scores from both sets. Application in ROC computation.

RELATED TO compute_ROCS

USAGE:
The script is designed to be called from the command line with the following arguments:
```bash
python ksm_scores_parser.py --pos <positive_scores> --neg <negative_scores> --out <output_file> --meth <method>
```

ARGUMENTS:
--pos, -p     : Path to the KSM scores file for positive sequences. Required.
--neg, -n     : Path to the KSM scores file for negative sequences. Required.
--out, -o     : Path to the output comparison file. Required.
--meth, -m    : Scoring method: "SUM", "BEST", or "MEAN". Required.

FUNCTIONS:
- averageOfList()    : Calculates the average of a list of numbers
- parse_ksm_scores() : Parses KSM score files and calculates summary statistics

INPUTS:
The input KSM score files should have the following format:
- Header line with "numSequence:<number>"
- Comment lines starting with "#"
- Data lines with tab-separated values including:
  - Motif information
  - Sequence ID
  - Score (in column 9, 0-based index 8)

OUTPUTS
The output file contains two columns:
- Positive sequence score
- Negative sequence score
Each row represents a matched pair of sequences (by rank) for the selected scoring method.

NOTES: 
- The script handles cases where sequences in the input files don't have scores by
  adding "fake" sequences with zero scores to ensure equal numbers of sequences.
- Three scoring methods are supported: SUM, BEST, and MEAN.
- The script assumes that both input files have the same number of sequences.

"""

# IMPORTS
from Bio import SeqIO
import sys
import os
import argparse
import operator

# LOAD ARGUMENTS
parser = argparse.ArgumentParser() 

parser.add_argument("--pos", "-p")
parser.add_argument("--neg", "-n")
parser.add_argument("--out", "-o")
parser.add_argument("--meth", "-m")
args = parser.parse_args()

## numSequence:1389
## numMotif:20
#Motif	SeqID	Motif_Name	SeqName	Match	SeqPos	Coord	Strand	Score

# FUNCTIONS: 
def averageOfList(num):
    sumOfNumbers = 0
    for t in num:
        sumOfNumbers = sumOfNumbers + t

    avg = sumOfNumbers / len(num)
    return avg

def parse_ksm_scores(FILE):
  f=open(FILE, "r")
  dico={}
  for l in f:
    if "numSequence" in l:
      numSeq=int(l.strip().split(":")[1])
      print("number of seq: %d\n" % numSeq)
      
    if l[0]!="#" and l[0:5]!="Motif":
      tmp=l.strip().split("\t")
      seq=tmp[3]
      score=float(tmp[8])
      if seq in dico.keys():
        dico[seq].append(score)
      else:
        dico[seq]=[score]  
  f.close()
  
  if len(dico)<numSeq:
    for i in range(numSeq-len(dico)):
      seq="fake%d" % i
      dico[seq]=[0]
  
  # sum score or take max or mean score
  dicoSum = {}
  dicoBest = {}
  dicoMean = {}
  for k in dico:
    dicoSum[k] = sum(dico[k])
    dicoBest[k] = max(dico[k])
    dicoMean[k] = averageOfList(dico[k])
  
  #sort from high to low score
  dicoSum_s = sorted(dicoSum.items(),key=operator.itemgetter(1),reverse=True)
  dicoBest_s = sorted(dicoBest.items(),key=operator.itemgetter(1),reverse=True)
  dicoMean_s = sorted(dicoMean.items(),key=operator.itemgetter(1),reverse=True)
  return {"SUM":dicoSum_s, "BEST":dicoBest_s, "MEAN":dicoMean_s}

# MAIN 
POS=parse_ksm_scores("%s" % args.pos)
NEG=parse_ksm_scores("%s" % args.neg)

# reminder: POS is dico of tuple 
# and POS["SUM"] is a tuple
size=min(len(POS["SUM"]), len(NEG["SUM"]))
print("a minimum of %d sequences had a kmer hit" % size)
cpt=0
fout=open("%s" % args.out, "w")
for i in range(size):
  scoPOS=POS[args.meth][i][1]
  scoNEG=NEG[args.meth][i][1]
  fout.write("%f\t%s\n" % (scoPOS, scoNEG))
fout.close()

