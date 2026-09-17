# -*- coding: utf-8 -*-
#!usr/bin/python

""" 
Motif Spacing Analyzer
======================
This script analyzes the spacing between transcription factor binding sites in DNA sequences. 
It identifies and categorizes different conformations of motif pairs (DR, ER, IR) and calculates 
the distances between them, considering various score thresholds.

RELATED TO compute_spacing2TFs

DESCRIPTION
1. Identifies three types of motif pair conformations:
	- DR (Direct Repeat): Both motifs on the same strand (-> -> or <- <-)
	- ER (Everted Repeat): Motifs on opposite strands, first motif on reverse strand (< ->)
	- IR (Inverted Repeat): Motifs on opposite strands, first motif on forward strand (-> <-)
2. Calculates distances between motif pairs within specified ranges.
3. Applies score thresholds to filter significant binding sites.
4. Handles both standard PWM and TFFM matrix outputs.
5. Generates detailed output files with spacing information.

INPUTS
- Positive file: A tab-separated file containing:
  * Header: Chromosome and coordinates (e.g., "chr1:1000-2000")
  * Position: Start position of the motif
  * Strand: "+", "-", "1", or "-1" indicating strand orientation
  * Score: Binding site score
  * Other optional columns
- Negative files: List of files containing sequences not bound by the transcription factor.

OUTPUTS
1. `<output_prefix>_spacing.tsv`: Contains information about motif pairs including:
2. `<output_prefix>_monomers.tsv`: Contains information about single motifs.
3. `<output_prefix>_counts.tsv`: Contains counts of motif pairs for each conformation and distance.
4. `<output_prefix>_enrichment.tsv`: Contains relative enrichment values for motif pair conformations.
5. `<output_prefix>_rate.tsv`: Contains absolute enrichment rates for motif pair conformations.
6. `<output_prefix>_NBofTFBS.txt`: Contains the number of transcription factor binding sites (TFBS) for each TF.

ARGUMENTS
- `--positive` or `-pos`: Path to the positive scores file.
- `--negative_sets` or `-neg`: List of paths to negative scores files.
- `--positive2` or `-pos2`: Path to the second positive scores file (optional).
- `--negative_sets2` or `-neg2`: List of paths to second negative scores files (optional).
- `--threshold` or `-th`: List of thresholds for the first TF (default: [-7, -8, -9, -10]).
- `--threshold2TF` or `-th2`: List of thresholds for the second TF (default: [0]).
- `--spacing_maxvalue` or `-smax`: Maximum spacing value (default: 30).
- `--spacing_minvalue` or `-smin`: Minimum spacing value (default: 0).
- `--offset_left` or `-ol`: Offset for the left side of the motif (default: 0).
- `--offset_right` or `-or`: Offset for the right side of the motif (default: 0).
- `--offset_left2TF` or `-ol2`: Offset for the left side of the second motif (default: 0).
- `--offset_right2TF` or `-or2`: Offset for the right side of the second motif (default: 0).
- `--len_motif` or `-lm`: Length of the first motif (default: 10).
- `--len_motif2TF` or `-lm2`: Length of the second motif (default: 10).
- `--output` or `-o`: Prefix for output files (default: "interdistances").
- `--write_inter` or `-wi`: Flag to write intermediate spacing results (default: False).
- `--no_absolute_panel` or `-nap`: Flag to disable absolute enrichment panel (default: False).

NOTES
- The script automatically detects whether TFFM matrices were used based on the threshold values.
- For TFFM matrices, the motif length is adjusted to account for the 2 bp extension needed for probability calculations.
- Dual-mode analysis is supported for analyzing two transcription factors simultaneously.
""" 

######################				General intels				######################				

#python_version  :2.7.9
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

def define_spacing(lib, spacing_max, spacing_min, offset_l, offset_r, len_motif, conformations, thresholds, write_inter, output,suffix):
	
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
	
	if write_inter:
		"""
		In case we need to report every spacing & conformation found.
		matricePosition1 & 2 are start position of the matrix for given score
		correctedPosition1 & 2 are start position of spacing distance (matrice position corrected by length of motif and possible offsets)
		"""
		with open(output+"_spacing_"+suffix+".tsv","w") as OUT:
			OUT.write("Peak" + "\t" + "Spacing" + "\t" + "Score1" + "\t" + "Score2" + "\t" + "matricePosition1" + "\t" + "matricePosition2" + "\t" + "correctedPosition1" + "\t" + "correctedPosition2" + "\n")
	minth=min(thresholds)
	for header in headers: # For each sequences
		sites = lib[np.where(lib['header'] == header)]
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
								if write_inter and th==minth: # report for every spacing/conformation found
									# report is done only for min threshold as you should be able to filter afterwards
									with open(output+"_spacing_"+suffix+".tsv","a") as OUT:
										OUT.write(header + "\t" + conformation + "_" + str(dist) + "\t" + str(elt1['score']) + "\t" + str(elt2['score']) + "\t" + str(elt1['position']) + "\t" + str(elt2['position']) + "\t" + str(position1) + "\t" + str(position2) + "\n")
	#for th in thresholds: # Debug report
		#print(str(th)+" - sum: "+str(conformations[th]['sum'])+"\nsum DR: "+str(conformations[th]['sumDR'])+"\nsumER: "+str(conformations[th]['sumER'])+"\nsumIR: "+str(conformations[th]['sumIR']))
	#for th in thresholds:
		#print("DR"+str(conformations[th]["DR"]))
		#print("ER"+str(conformations[th]["ER"]))
		#print("IR"+str(conformations[th]["IR"]))
	return conformations

def  define_spacing_2TF(lib, lib_2, spacing_max, spacing_min, offset_l, offset_r, offset_l2, offset_r2, len_motif, len_motif2,conformations, thresholds, thresholds_2, write_inter, output,suffix):
	
	headers = np.unique(np.sort(lib['header'])) # Getting all sequences to analyses.
	
	# preparing dictionnary to contain spacing results
	if len(conformations.keys()) == 0:
		thresholds_ref=thresholds
		for th in thresholds_ref:
			conformations[th]={}
			conformations[th]['ER']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['IR']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['DR12']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['DR21']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['sum']=0
			conformations[th]['sumER']=0
			conformations[th]['sumIR']=0
			conformations[th]['sumDR12']=0
			conformations[th]['sumDR21']=0
	if write_inter:
		"""
		In case we need to report every spacing & conformation found.
		matricePosition1 & 2 are start position of the matrix for given score
		correctedPosition1 & 2 are start position of spacing distance (matrice position corrected by length of motif and possible offsets)
		"""
		with open(output+"_spacing_"+suffix+".tsv","w") as OUT:
			OUT.write("Peak" + "\t" + "Spacing" + "\t" + "Score1" + "\t" + "Score2" + "\t" + "matricePosition1" + "\t" + "matricePosition2" + "\t" + "correctedPosition1" + "\t" + "correctedPosition2" + "\n")
	minth=min(thresholds)
	minth_2=min(thresholds_2)
	for header in headers: # For each sequences
		sites = lib[np.where(lib['header'] == header)]
		sites_2 = lib_2[np.where(lib_2['header'] == header)]
		sites = np.sort(sites,order='position')
		sites_2 = np.sort(sites_2,order='position')
		if len(sites) > 0 and len(sites_2) > 0 : # if sequences contain at least 2 sites
			for x1 in range(0,len(sites)):
				elt1=sites[x1]
				for x2 in range(x1,len(sites_2)):
					elt2=sites_2[x2]
					conformation="NA"
					dist = "NA"
					if elt1["strand"] == elt2["strand"]: # DR
						if elt1["strand"] == "1" or elt1["strand"] == "+": #DR 1> 2>
							dist = (elt2["position"] + offset_l2) - (elt1["position"] + len_motif - offset_r)
							position1=elt1["position"] - 1
							position2=elt2["position"] + len_motif2
							if "12" in suffix:
								conformation = "DR12"
							elif "21" in suffix:
								conformation = "DR21"
						else: #DR <1 <2
							dist = (elt2["position"] + offset_r2) - (elt1["position"] + len_motif - offset_l)
							position1=elt1["position"] - 1
							position2=elt2["position"] + len_motif2
							if "12" in suffix:
								conformation = "DR21"
							elif "21" in suffix:
								conformation = "DR12"
						#conformation = "DR"
					if "12" in suffix: # avoiding to register same ER and IR twice as they are the same between suffix 12 and 12
						if (elt1["strand"] == "-1" and elt2["strand"] =="1") or (elt1["strand"] == "-" and elt2["strand"] =="+"): #ER <1 2>
							dist = (elt2["position"] + offset_l2) - (elt1["position"] + len_motif - offset_l)
							conformation = "ER"
							position1=elt1["position"] - 1
							position2=elt2["position"] + len_motif2
						if (elt1["strand"] == "1" and elt2["strand"] =="-1") or (elt1["strand"] == "+" and elt2["strand"] =="-"): #IR 1> <2
							dist = (elt2["position"] + offset_r2) - (elt1["position"] + len_motif - offset_r)
							conformation = "IR"
							position1=elt1["position"] - 1
							position2=elt2["position"] + len_motif2
					if dist!= "NA" and (spacing_min <= dist <= spacing_max):
						for th,th2 in zip(thresholds, thresholds_2):
							if elt1['score'] > th and elt2['score'] > th2: # making sure each sites have score higher than thresholds and add to corresponding sums
								if suffix=="pos12":
									conformations[th][conformation][dist]+=1
									conformations[th]['sum']+=1
									conformations[th]['sum'+str(conformation)]+=1
								else:
									conformations[th2][conformation][dist]+=1
									conformations[th2]['sum']+=1
									conformations[th2]['sum'+str(conformation)]+=1
								if write_inter and th==minth: # report for every spacing/conformation found
									# report is done only for min threshold as you should be able to filter afterwards
									
									with open(output+"_spacing_"+suffix+".tsv","a") as OUT:
										OUT.write(header + "\t" + conformation + "_" + str(dist) + "\t" + str(elt1['score']) + "\t" + str(elt2['score']) + "\t" + str(elt1['position']) + "\t" + str(elt2['position']) + "\t" + str(position1) + "\t" + str(position2) + "\n")
	#for th in thresholds: # Debug report
		#print(str(th)+" - sum: "+str(conformations[th]['sum'])+"\nsum DR: "+str(conformations[th]['sumDR'])+"\nsumER: "+str(conformations[th]['sumER'])+"\nsumIR: "+str(conformations[th]['sumIR']))
	#for th in thresholds:
		#print("DR"+str(conformations[th]["DR"]))
		#print("ER"+str(conformations[th]["ER"]))
		#print("IR"+str(conformations[th]["IR"]))
	
	if suffix=="pos12":
		with open(output+"_counts.tsv","w") as OUT2:
			for th in thresholds:
				OUT2.write(str(th)+"\n")
				for conform in conformations[th].keys():
					if "sum" not in conform:
						for distance in range(len(conformations[th][conform])):
							OUT2.write(str(conform)+"\t"+str(distance)+"\t"+str(conformations[th]["DR21"][distance])+"\n")
	else: # "pos21"
		with open(output+"_counts.tsv","a") as OUT2:
			for th in thresholds_2:
				OUT2.write(str(th)+"\n")
				for distance in range(len(conformations[th]["DR21"])):
					OUT2.write(str("DR21")+"\t"+str(distance)+"\t"+str(conformations[th]["DR21"][distance])+"\n")
	return conformations
						
############################
def write_enrichment(relative_conformations, spacing_max, spacing_min, thresholds, output, no_absolute_panel, rate):
	dist_list = np.arange(spacing_min,spacing_max + 1)
	with open(output+"_enrichment.tsv","w") as OUT:
		OUT.write("threshold\tconformation\tspace\tenrichment\n")
		for th in thresholds:
			for conformation in ['DR', 'IR', 'ER']:
				for dist in dist_list:
					to_print=True
					if relative_conformations['neg'][th][conformation][dist]==0:
						print("======\nthreshold: "+str(th)+"\nconformation: "+str(conformation)+"\ndistance: "+str(dist)+"\nNo sites found in this conformation for negative set\nThreshold set too high\n=====")
						to_print=False
					if relative_conformations['pos'][th][conformation][dist]==0:
						print("======\nthreshold: "+str(th)+"\nconformation: "+str(conformation)+"\ndistance: "+str(dist)+"\nNo sites found in this conformation for positive set\nThreshold set too high\n=====")
						to_print=False
					if to_print:
						OUT.write(str(th) + "\t" + conformation + "\t" + str(dist) + "\t" + str(relative_conformations['pos'][th][conformation][dist]/relative_conformations['neg'][th][conformation][dist])+"\n")
			for dist in dist_list:
				if relative_conformations['neg'][th]['ALL'][dist] != 0:
					OUT.write(str(th) + "\t" + 'all' + "\t" + str(dist) + "\t" + str(relative_conformations['pos'][th]['ALL'][dist]/relative_conformations['neg'][th]['ALL'][dist])+"\n")
				else:
					OUT.write(str(th) + "\t" + 'all' + "\t" + str(dist) + "\t" + "0.0" +"\n")
	if no_absolute_panel:
		pass
	else:
		with open(output+"_rate.tsv","w") as OUT:
			OUT.write("threshold\tconformation\trate\n")
			i=0
			for th in thresholds:
				for conformation,indices in zip(['DR', 'IR', 'ER'],[0,2,1]):
					OUT.write(str(th) + "\t" + conformation + "\t" + str(rate[i][indices])+"\n")
				i+=1


######################				MAIN				######################
def main(	positive_file,
			negative_files,
			thresholds,
			spacing_max,
			spacing_min,
			offset_l,
			offset_r,
			len_motif,
			output,
			write_inter,
			no_absolute_panel,
			positive_file_2,
			negative_files_2,
			thresholds_2,
			offset_l2,
			offset_r2,
			len_motif2
		):
	th_min = min(thresholds)
	th_min_2 = min(thresholds_2)
	#print(thresholds) # Debug print
	tffm=False
	dualmode=False
	tffm2=False
	# if at least one th is greater than 0, then tffm matrix have been used
	for th in thresholds:
		if th >0:
			tffm=True
	if positive_file_2[0] != "NA":
		dualmode=True
		if len(negative_files) != len(negative_files_2):
			print("ERROR, number of negative files used for TF1 and TF2 are different")
			sys.exit(1)
		if th in thresholds_2:
			if th >0:
				tffm2=True
	
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
		#positive_set = np.genfromtxt(positive_file, dtype=[('header','S300'),('position',int),('strand','S5'),('score',float)])
		length_positive=compute_Length_seq(positive_set, len_motif)
	
	
	if tffm2 and dualmode:
		positive_set_2 = np.genfromtxt(positive_file_2, dtype=[('header','U300'), ('position',int), ('stop',int), ('strand','U5'), ('seq','U300'), ('type','U10'),('lengthM',int), ('score',float)])
		len_motif2 = (positive_set_2['stop'][1] - positive_set_2['position'][1] + 1)
		
		length_positive_2=compute_Length_seq(positive_set, len_motif2+2)
	elif not tffm2 and dualmode:
		#positive_set_2 = np.genfromtxt(positive_file_2, dtype=[('header','S300'),('position',int),('strand','S5'),('score',float)])
		positive_set_2 = np.genfromtxt(positive_file_2, dtype=[('header','U300'), ('position',int), ('stop',int), ('strand','U5'), ('seq','U300'), ('type','U10'),('lengthM',int), ('score',float)])
		length_positive_2=compute_Length_seq(positive_set_2, len_motif2)
		
	# we filter out sequences & sites with score < minimum threhold to avoid unecessary computation (gain of time)
	#print("number of sequences to analyse before filtration by threshold: " + str(len(np.unique(np.sort(positive_set['header'])))))
	positive_set = positive_set[np.where( positive_set['score'] > th_min )]
	if dualmode:
		positive_set_2 = positive_set_2[np.where( positive_set_2['score'] > th_min_2 )]
	#print("number of sequences to analyse after filtration by threshold: " + str(len(np.unique(np.sort(positive_set['header'])))))
	with open(output+"_NBofTFBS.txt","w") as OUT:
		OUT.write("TFBS_for_TF1\t"+str(len(positive_set['header']))+"\n")
		OUT.write("TFBS_for_TF2\t"+str(len(positive_set_2['header']))+"\n")
	conformations={}
	
	# finding ER, DR and IR conformation for each spacing
	print("defining spacing for positive set")
	conformations['pos']={}
	if not dualmode:
		conformations['pos'] = define_spacing(positive_set, spacing_max, spacing_min, offset_l, offset_r, len_motif,conformations['pos'], thresholds, write_inter, output,"pos")
	else:
		conformations['pos'] = define_spacing_2TF(positive_set, positive_set_2, spacing_max, spacing_min, offset_l, offset_r, offset_l2, offset_r2, len_motif, len_motif2,conformations['pos'], thresholds, thresholds_2, write_inter, output,"pos12")
		
		conformations['pos'] = define_spacing_2TF(positive_set_2, positive_set, spacing_max, spacing_min, offset_l2, offset_r2, offset_l, offset_r, len_motif2, len_motif,conformations['pos'], thresholds_2, thresholds, write_inter, output,"pos21")
	
	if not write_inter:
		conformations['neg']={}
		i=0
		length_negative=0
		length_negative_2=0
		if not dualmode:
			negative_files_2=negative_files
		for elt,elt2 in zip(negative_files, negative_files_2):
			i+=1
			print("defining spacing for negative set "+str(i))
			if tffm:
				negative_set = np.genfromtxt(elt, dtype=[('header','S300'), ('position',int), ('stop',int), ('strand','S5'), ('seq','S300'), ('type','S10'),('lengthM',int), ('score',float)])
				"""
				length of motif is extended for maximum number of sites computation. See explanation done for positive set )
				"""
				length_negative+=compute_Length_seq(negative_set, len_motif+2)
			else:
				negative_set = np.genfromtxt(elt, dtype=[('header','S300'),('position',int),('strand','S5'),('score',float)])
				length_negative+=compute_Length_seq(negative_set, len_motif)
			if tffm2 and dualmode:
				negative_set_2 = np.genfromtxt(elt2, dtype=[('header','S300'), ('position',int), ('stop',int), ('strand','S5'), ('seq','S300'), ('type','S10'),('lengthM',int), ('score',float)])
				length_negative_2+=compute_Length_seq(negative_set_2, len_motif2+2)
			elif dualmode:
				negative_set_2 = np.genfromtxt(elt2, dtype=[('header','S300'),('position',int),('strand','S5'),('score',float)])
				length_negative_2+=compute_Length_seq(negative_set_2, len_motif2)
			negative_set = negative_set[np.where( negative_set['score'] > th_min )]
			if dualmode:
				negative_set_2 = negative_set_2[np.where( negative_set_2['score'] > th_min_2 )]
			#print("number of sequences to analyse after filtration by threshold: " + str(len(np.unique(np.sort(negative_set['header'])))))
			# finding ER, DR and IR conformation for each spacing & saving in one dictionnary for all negative set.
			if not dualmode:
				conformations['neg'] = define_spacing(negative_set, spacing_max, spacing_min, offset_l, offset_r, len_motif, conformations['neg'], thresholds, write_inter, output,"neg"+str(i))
			elif dualmode:
				conformations['neg'] = define_spacing_2TF(negative_set, negative_set_2, spacing_max, spacing_min, offset_l, offset_r, offset_l2, offset_r2, len_motif, len_motif2, conformations['neg'], thresholds, thresholds_2, write_inter, output,"neg"+str(i)+"12")
				
				conformations['neg'] = define_spacing_2TF(negative_set_2, negative_set, spacing_max, spacing_min, offset_l2, offset_r2, offset_l, offset_r, len_motif2, len_motif, conformations['neg'], thresholds_2, thresholds, write_inter, output,"neg"+str(i)+"21")
	
		# normalisation computation to take into account the higher number of sequences treated for the negative set
		#print(length_positive)
		#print(length_negative)
		if not dualmode:
			norm = length_positive/float(length_negative)
		else: # dualmode
			norm = max([length_positive/float(length_negative),length_positive_2/float(length_negative_2)])
		print("normalisation: "+str(norm))
		#print(conformations) # Debug print
		
		rate=[]
		for th in thresholds:
			"""
			Absolute Enrichment computation:
			( sum(DRpos) / sum(DRneg) ) * all_possible_site_pos/all_possible_site_neg
			"""
			rate.append(	[divide(conformations['pos'][th]['sumDR'],conformations['neg'][th]['sumDR'])/norm,
							divide(conformations['pos'][th]['sumER'],conformations['neg'][th]['sumER'])/norm,
							divide(conformations['pos'][th]['sumIR'],conformations['neg'][th]['sumIR'])/norm ])
		#print(rate) # Debug print
		
		relative_conformations={}
		relative_conformations['pos']={}
		relative_conformations['neg']={}
		"""
		Relative Enrichment computation: 
		sum(DRpos th|dist) / ( sum(NRpos) ) / sum(DRneg th|dist) / ( sum(NRneg) )
		
		NR : ER + DR + IR
		
		"""
		for key in relative_conformations.keys():
			for th in thresholds:
				relative_conformations[key][th]={}
				relative_conformations[key][th]['ER']=[ divide(x, float(conformations[key][th]['sum'])) for x in conformations[key][th]['ER'] ]
				relative_conformations[key][th]['IR']=[ divide(x, float(conformations[key][th]['sum'])) for x in conformations[key][th]['IR'] ]
				relative_conformations[key][th]['DR']=[ divide(x, float(conformations[key][th]['sum'])) for x in conformations[key][th]['DR'] ]
				relative_conformations[key][th]['ALL']=[ divide(x+y+z, float(conformations[key][th]['sum'])) for x,y,z in zip(conformations[key][th]['DR'],conformations[key][th]['ER'],conformations[key][th]['IR']) ]
		#for th in thresholds: # Debug prints
			# for conformation in ['DR', 'ER','IR','ALL']
				#print(conformation+" pos relative enrichment "+str(relative_conformations['pos'][th][conformation]))
				#print(conformation+" neg relative enrichment "+str(relative_conformations['neg'][th][conformation]))
		write_enrichment(relative_conformations, spacing_max, spacing_min, thresholds, output, no_absolute_panel, rate)
	
######################				Parsing arguments				######################
parser = argparse.ArgumentParser()
parser.add_argument("--positive", "-pos", help='scores file of sequences bound by TF, generated by scores.py')
parser.add_argument("--negative_sets", "-neg", nargs='+', help='list of scores files of sequences not bound by TF, generated by scores.py' )
parser.add_argument("--positive2", "-pos2", help='scores file of sequences bound by TF, generated by scores.py', default=["NA"])
parser.add_argument("--negative_sets2", "-neg2", nargs='+', help='list of scores files of sequences not bound by TF, generated by scores.py', default="NA")
parser.add_argument("--threshold", "-th",nargs='+',type = float, default= [-7, -8, -9, -10])
parser.add_argument("--threshold2TF", "-th2",nargs='+',type = float, default= [0])
parser.add_argument("--spacing_maxvalue", "-smax", type = int, default= 30)
parser.add_argument("--spacing_minvalue", "-smin", type = int, default= 0)
parser.add_argument("--offset_left", "-ol", type = int, default=0)
parser.add_argument("--offset_right", "-or", type = int, default=0)
parser.add_argument("--offset_left2TF", "-ol2", type = int, default=0)
parser.add_argument("--offset_right2TF", "-or2", type = int, default=0)
parser.add_argument("--len_motif", "-lm", type = int, default=10)
parser.add_argument("--len_motif2TF", "-lm2", type = int, default=10)
parser.add_argument("--output", "-o", type = str, default="interdistances")
parser.add_argument("--write_inter", "-wi",action='store_true', default= False)
parser.add_argument("--no_absolute_panel", "-nap",action='store_true', default= False)


args = parser.parse_args()

#print(args.negative_sets2) #Debug print
# print(args.offset_right2TF)
# print(args.offset_right)
main(	args.positive,
		args.negative_sets,
		args.threshold,
		args.spacing_maxvalue,
		args.spacing_minvalue,
		args.offset_left,
		args.offset_right,
		args.len_motif,
		args.output,
		args.write_inter,
		args.no_absolute_panel,
		args.positive2,
		args.negative_sets2,
		args.threshold2TF,
		args.offset_left2TF,
		args.offset_right2TF,
		args.len_motif2TF
	)
	
#python get_interdistances.py -pos /home/304.6-RDF/Jeremy/Arf_Paper_fig1_maker/results/2018_03_29_ARF2_fullset/ARF2_MP/Find_motifs/ARF2/scores/testing_pos_set.fasta.scores -neg /home/304.6-RDF/Jeremy/Arf_Paper_fig1_maker/results/2018_03_29_ARF2_fullset/ARF2_MP/Find_motifs/ARF2/scores/testing_neg_set.fasta.scores -ol 2 -or 2


