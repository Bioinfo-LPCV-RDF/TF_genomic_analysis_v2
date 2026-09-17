#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Probability of Occupancy (POcc) Calculator
==========================================

DECSRIPTION
This script calculates the Probability of Occupancy (POcc) for transcription factor
binding sites based on scoring data. It implements the formula from Moyroud et al. (2011)
to convert raw binding scores into probabilities of occupancy.

The script reads a file containing scores for each binding site, calculates the POcc
for each gene by summing the probabilities of its individual binding sites, and writes
the results to an output file.

RELATED TO compute_ROCS

USAGE
The script is designed to be called from the command line with the following arguments:
```bash
python pocc_calculator.py --scores <scores_file> --output <output_file> [--a <a_value>] [--b <b_value>]
``` 

ARGUMENTS:
--scores, -s   : Path to the scores file. Required.
--output, -o   : Path to the output file. Required.
--a, -a        : Parameter 'a' for the POcc formula. Default: 1.0.
--b, -b        : Parameter 'b' for the POcc formula. Default: 0.0.

FUNCTIONS: 
Ka(score, a, b) : Converts a raw score to a binding probability using the formula: exp(-(b-score)/a)

INPUT 
The input scores file should be a tab-separated file with at least 8 columns:
- Column 1: Gene identifier
- Column 8: Binding site score (0-based index 7)

Example: 
Gene1   chr1   100   110   ...   0.85
Gene1   chr1   150   160   ...   0.92
Gene2   chr2   200   210   ...   0.78

OUTPUT:
The output file contains one POcc value per line, each corresponding to a gene from the input file.

Example Output:
1.2345
0.8765
...

NOTES: 
- The default parameters (a=1.0, b=0.0) are placeholders. For real analysis, these should
  be set to values appropriate for your specific transcription factor (e.g., a=2.5663,
  b=0.3598 for LFY as in Moyroud et al., 2011).
- The script assumes the input file is sorted by gene.
- Each gene's POcc is calculated as the sum of probabilities for all its binding sites.
- The script handles multiple binding sites per gene by accumulating their probabilities.

LEAFY CONSTANS: (according to Moyroud and al., 2011 - valid for SYM, full site matrix LFY.pfm)
- a=2.5663
- b=0.3598
- TF = exp(b/a)
- threshold = -20
- output = "/home/304.6-RDF/Pauline/POcc" 
- matrix = "/home/304.6-RDF/data/LFY.pfm"
- sequences = "/home/304.6-RDF/FParcy_prom_results/apetala3/promoters_all_clade/promoters_seq_by_length/prom1000.fa"

""" 

# IMPORTS
from math import exp
import os
import argparse

# FUNCTIONS
def Ka(score,a,b) :
	return exp(-(b-float(score))/a)

# MAIN ##################
def main(scores,output,a,b,writegene) :
	
	a = float(a)
	b = float(b)

	TF = exp(b/a)
	
	with open(scores, "r") as score_file, open(output,"w") as pocc_scores :
		s = 0
		current_gene = "" 
		for line in score_file :
			line = line.split('\t')
			# print(line[0])
			if current_gene != "" and line[0] == current_gene : 
				s+= Ka(line[7],a,b)*TF/(1+Ka(line[7],a,b)*TF) 
			else :
				if writegene:
					pocc_scores.write(current_gene+str(s)+"\n")
				else :
					pocc_scores.write(str(s)+"\n") #on est passés à une autre espèce donc on écrit la POcc pour la précedente
				s = Ka(line[7],a,b)*TF/(1+Ka(line[7],a,b)*TF) #on recommence à calculer la POcc pour la nouvelle espèce
				current_gene = line[0]
		if writegene:
			pocc_scores.write(current_gene+str(s)+"\n")
		else :
			pocc_scores.write(str(s)+"\n") 
		
	print('[INFO] - POcc computation done')
# ########################

# ARGUMENTS
if __name__ == "__main__": 
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--scores', help = 'scores files produced by ', required = True)
	parser.add_argument('-o', '--output', help = 'output path', required = True)
	parser.add_argument('-a', '--a', help = '', default=1.0)
	parser.add_argument('-b', '--b', help = '', default=0.0)
	parser.add_argument("--writegene", "-w",action='store_true', default= False)
	args = parser.parse_args()
	main(args.scores, args.output,args.a,args.b,args.writegene)



