#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Motif Spacing Analyzer
======================

DESCRIPTION:
This Python script analyzes the spacing between transcription factor binding sites
in DNA sequences. It identifies and categorizes different conformations of motif pairs
(DR, ER, IR) and calculates the distances between them, considering various score
thresholds.

RELATED TO compute_space

FEATURES:
- Identifies three types of motif pair conformations:
  * DR (Direct Repeat): Both motifs on the same strand (-> -> or <- <-)
  * ER (Everted Repeat): Motifs on opposite strands, first motif on reverse strand (< ->)
  * IR (Inverted Repeat): Motifs on opposite strands, first motif on forward strand (-> <-)

- Calculates distances between motif pairs within specified ranges
- Applies score thresholds to filter significant binding sites
- Handles both standard PWM and TFFM matrix outputs
- Generates detailed output files with spacing information

USAGE:
The script is designed to be called from the command line with the following arguments:
```bash
python motif_spacing.py --positive <positive_file> [--threshold <thresholds>]
                        [--spacing_maxvalue <max_spacing>] [--spacing_minvalue <min_spacing>]
                        [--offset_left <left_offset>] [--offset_right <right_offset>]
                        [--len_motif <motif_length>] [--output <output_prefix>]
``` 

ARGUMENTS:
--positive, -pos          : Scores file of sequences bound by TF. Required
--threshold, -th          : Score thresholds (multiple values allowed). Default: [-7, -8, -9, -10]
--spacing_maxvalue, -smax : Maximum spacing between motifs. Default: 30
--spacing_minvalue, -smin : Minimum spacing between motifs. Default: 0
--offset_left, -ol        : Left offset for motif positioning. Default: 0
--offset_right, -or       : Right offset for motif positioning. Default: 0
--len_motif, -lm          : Length of the motif. Default: 10
--output, -o              : Output file prefix. Default: "interdistances"

FUNCTIONS: 
- divide(a, b)                       : Safe division function that returns 0 if division by zero occurs
- compute_Length_seq(lib, len_motif) : Calculates the total length of sequences to analyze
- define_spacing()                   : Main function that identifies motif pairs and calculates spacings

INPUT FILE: 
The positive file should be a tab-separated file containing:
- Header     : Chromosome and coordinates (e.g., "chr1:1000-2000")
- Position   : Start position of the motif
- Strand     : "+", "-", "1", or "-1" indicating strand orientation
- Score      : Binding site score
- Other optional columns

OUTPUTS: 
- <output_prefix>_spacing.tsv: Contains information about motif pairs including:
 - Chromosome coordinates
 - Conformation type (DR, ER, IR)
 - Spacing distance
 - Individual motif scores
 - Motif positions
- <output_prefix>_monomers.tsv: Contains information about single motifs

NOTES: 
- The script automatically detects whether TFFM matrices were used based on the
  threshold values (if any threshold > 0)
- For TFFM matrices, the motif length is automatically adjusted to account for
  the 2 bp extension needed for probability calculations
"""

######################				import packages				######################				

import argparse #v1.1
import numpy as np # 1.16.1
import sys
import math
import re # 2.2.1

######################				Functions				######################
def divide(a, b):
	try:
		return a/float(b)
	except:
		return 0.0

def compute_Length_seq(lib,len_motif):
	headers = np.unique(np.sort(lib['header']))
	length_set=0
	for header in headers:
		lengthseq=int(header.split(":")[-1].split("-")[1])-int(header.split(":")[-1].split("-")[0]) + 1 - len_motif
		length_set+=lengthseq
	return length_set


def define_spacing(lib, spacing_max, spacing_min, offset_l, offset_r, len_motif, conformations, thresholds, output):
	
	headers = np.unique(np.sort(lib['header'])) # Getting all sequences to analyses.
	
	# preparing dictionnary to contain spacing results
	if len(conformations.keys()) == 0:
		for th in thresholds:
			conformations[th]={}
			conformations[th]['ER']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['IR']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['DR']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['sum']=0
			conformations[th]['sumER']=0
			conformations[th]['sumIR']=0
			conformations[th]['sumDR']=0
		"""
		In case we need to report every spacing & conformation found.
		matricePosition1 & 2 are start position of the matrix for given score
		correctedPosition1 & 2 are start position of spacing distance (matrice position corrected by length of motif and possible offsets)
		"""
	
			
	lines_to_write=[]
	monomer=[]
	#used_header=[]
	minth=min(thresholds)

	for header in headers: # For each sequences
		sites = lib[np.where(lib['header'] == header)]
		chrom=header.split(":")[0]
		start=header.split(":")[1].split("-")[0]
		stop=header.split(":")[1].split("-")[1]
		header=chrom+"\t"+start+"\t"+stop
		sites = np.sort(sites,order='position')
		if len(sites) > 1: # if sequences contain at least 2 sites
			for x1 in range(0,len(sites)-1):
				elt1=sites[x1]
				for x2 in range(x1+1,len(sites)):
					elt2=sites[x2]
					"""
					double for loop to analyses possibles combinaisons of Spacing.
					DR: define by same strandness of both sites
					ER: first sites' strand must be forward and second sites' strand must be reverse
					IR: first sites' strand must be reverse and second sites' strand must be forward
					
					note that distance computation is dependant on conformation (offset used are not the same)
					
					"""
					conformation="NA"
					dist = "NA"
					if elt1["strand"] == elt2["strand"]: # DR
						if elt1["strand"] == "1" or elt1["strand"] == "+": #DR > >
							dist = (elt2["position"] + offset_l) - (elt1["position"] + len_motif - offset_r)
							position1=elt1["position"] + len_motif - offset_r
							position2=elt2["position"] + offset_l
						else: #DR < <
							dist = (elt2["position"] + offset_r) - (elt1["position"] + len_motif - offset_l)
							position1=elt1["position"] + len_motif - offset_l
							position2=elt2["position"] + offset_r
						conformation = "DR"
					if (elt1["strand"] == "-1" and elt2["strand"] =="1") or (elt1["strand"] == "-" and elt2["strand"] =="+"): #ER < >
						dist = (elt2["position"] + offset_l) - (elt1["position"] + len_motif - offset_l)
						conformation = "ER"
						position1=elt1["position"] + len_motif - offset_l
						position2=elt2["position"] + offset_l
					if (elt1["strand"] == "1" and elt2["strand"] =="-1") or (elt1["strand"] == "+" and elt2["strand"] =="-"): #IR > <
						dist = (elt2["position"] + offset_r) - (elt1["position"] + len_motif - offset_r)
						conformation = "IR"
						position1=elt1["position"] + len_motif - offset_r
						position2=elt2["position"] + offset_r
					if dist!= "NA" and (spacing_min <= dist <= spacing_max):
						for th in thresholds:
							if elt1['score'] > th and elt2['score'] > th: # making sure each sites have score higher than thresholds and add to corresponding sums
								conformations[th][conformation][dist]+=1
								conformations[th]['sum']+=1
								conformations[th]['sum'+str(conformation)]+=1
								if th==minth: # report for every spacing/conformation found
									# report is done only for min threshold as you should be able to filter afterwards
									lines_to_write.append(header + "\t" + conformation + "\t" + str(dist) + "\t" + str(elt1['score']) + "\t" + str(elt2['score']) + "\t" + str(elt1['position']) + "\t" + str(elt2['position']) + "\t" + str(position1) + "\t" + str(position2))
									#used_header.append(header)
					if elt1['score'] > th and dist=="NA":
						monomer.append(header+"\t"+"monomer"+"\t"+"0"+"\t"+str(elt1['score'])+"\t"+str(elt1['position']))
		#if header not in used_header:
			#lines_to_write.append(header + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")
	#for th in thresholds: # Debug report
		#print(str(th)+" - sum: "+str(conformations[th]['sum'])+"\nsum DR: "+str(conformations[th]['sumDR'])+"\nsumER: "+str(conformations[th]['sumER'])+"\nsumIR: "+str(conformations[th]['sumIR']))
	#for th in thresholds:
		#print("DR"+str(conformations[th]["DR"]))
		#print("ER"+str(conformations[th]["ER"]))
		#print("IR"+str(conformations[th]["IR"]))

	with open(output+"_spacing.tsv","w") as OUT:
		OUT.write("chr"+ "\t" + "Start" + "\t" + "End" + "\t" + "Conf" + "\t" + "Space" + "\t" + "Score1" + "\t" + "Score2" + "\t" + "matricePosition1" + "\t" + "matricePosition2" + "\t" + "correctedPosition1" + "\t" + "correctedPosition2" + "\n")
		OUT.write("\n".join(lines_to_write))
	with open(output+"_monomers.tsv","w") as OUT2:
		OUT2.write("chr"+ "\t" + "Start" + "\t" + "End" + "\t" + "Conf" + "\t" + "Space" + "\t" + "Score1" + "\t" + "matricePosition1" + "\n")
		OUT2.write("\n".join(monomer))
	return conformations

######################				MAIN				######################
def main(	positive_file,
			thresholds,
			spacing_max,
			spacing_min,
			offset_l,
			offset_r,
			len_motif,
			output,
		):
	th_min = min(thresholds)
	#print(thresholds) # Debug print
	tffm=False
	# if at least one th is greater than 0, then tffm matrix have been used
	for th in thresholds:
		if th > 0:
			tffm=True
	if tffm:
		positive_set = np.genfromtxt(positive_file, dtype=[('header','U300'), ('position',int), ('stop',int), ('strand','U5'), ('seq','U300'), ('type','U10'),('lengthM',int), ('score',float)])
		len_motif = (positive_set['stop'][1] - positive_set['position'][1] + 1)
		"""
		length of motif is extended for maximum number of sites computation because TFFM matrices needs 2 bp in front of the matrix to take in consideration probability of apparition of the site.
		We don't need this correction for spacing computation since score are reported with correct length of motif (length similar to PWM's length)
		"""
		length_positive=compute_Length_seq(positive_set, len_motif+2)
	else:
		positive_set = np.genfromtxt(positive_file, dtype=[('header','U300'), ('position',int), ('stop',int), ('strand','U5'), ('seq','U300'), ('type','U10'),('lengthM',int), ('score',float)])
		length_positive=compute_Length_seq(positive_set, len_motif)
	
	# we filter out sequences & sites with score < minimum threhold to avoid unecessary computation (gain of time)
	#print("number of sequences to analyse before filtration by threshold: " + str(len(np.unique(np.sort(positive_set['header'])))))
	positive_set = positive_set[np.where( positive_set['score'] > th_min )]
	#print("number of sequences to analyse after filtration by threshold: " + str(len(np.unique(np.sort(positive_set['header'])))))
	conformations={}
	
	# finding ER, DR and IR conformation for each spacing
	print("defining spacing for positive set")
	conformations['pos']={}
	conformations['pos'] = define_spacing(positive_set, spacing_max, spacing_min, offset_l, offset_r, len_motif,conformations['pos'], thresholds, output)

######################				Parsing arguments				######################

parser = argparse.ArgumentParser()
parser.add_argument("--positive", "-pos", help='scores file of seqences bound by TF, generated by scores.py')
parser.add_argument("--threshold", "-th",nargs='+',type = float, default= [-7, -8, -9, -10])
parser.add_argument("--spacing_maxvalue", "-smax", type = int, default= 30)
parser.add_argument("--spacing_minvalue", "-smin", type = int, default= 0)
parser.add_argument("--offset_left", "-ol", type = int, default=0)
parser.add_argument("--offset_right", "-or", type = int, default=0)
parser.add_argument("--len_motif", "-lm", type = int, default=10)
parser.add_argument("--output", "-o", type = str, default="interdistances")

args = parser.parse_args()

#print(args.negative_sets) #Debug print
main(	args.positive,
		args.threshold,
		args.spacing_maxvalue,
		args.spacing_minvalue,
		args.offset_left,
		args.offset_right,
		args.len_motif,
		args.output,
	)
	
#python get_interdistances.py -pos /home/304.6-RDF/Jeremy/Arf_Paper_fig1_maker/results/2018_03_29_ARF2_fullset/ARF2_MP/Find_motifs/ARF2/scores/testing_pos_set.fasta.scores -neg /home/304.6-RDF/Jeremy/Arf_Paper_fig1_maker/results/2018_03_29_ARF2_fullset/ARF2_MP/Find_motifs/ARF2/scores/testing_neg_set.fasta.scores -ol 2 -or 2


