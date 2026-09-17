# -*- coding: utf-8 -*-
# !usr/bin/python

################################################			   General intels				################################################			   
# description	 : 
# author		  : 
# date			: 
# version		 :
# usage		   : 
# python_version  : 2.7.9
################################################			   import packages				################################################			   
import re
import tempfile 
import os
import argparse # version 1.1
import pybedtools # version 0.7.10
from Bio import SeqIO # version 1.70
import numpy as np # version 1.8.2
import matplotlib.pyplot as plt
import sklearn.metrics
from string import maketrans
import sys
################################################			   Parsing arguments				################################################
parser = argparse.ArgumentParser(description = "") 
parser.add_argument("--data","-d", help = "positive set first, negative set after",nargs = 2, required = True)
parser.add_argument("--genome","-g", help = "", required = True)
parser.add_argument("--matrix","-m", nargs = "*", default = [], help = "matrix to use")
parser.add_argument("--threhold","-t", nargs = "*", help = "")
########### saving arguments to variables
args = parser.parse_args()
(posSet, negSet) = args.data
matrix = args.matrix
threhold = args.threhold
genome = args.genome
########### Naming temporary files
temp_var = tempfile.NamedTemporaryFile(dir = ".")
pos_temp_filename = temp_var.name
temp_var.close()

temp_var = tempfile.NamedTemporaryFile(dir = ".")
neg_temp_filename = temp_var.name
temp_var.close()

################################################			   Functions				################################################
"""
from a bedfile, creates Fasta object containing all the sequences in file
"""
def make_fasta(filename, tempfile, genome):
	bedfile = pybedtools.BedTool(filename).saveas(tempfile) # reading file and saving it in a temp file
	bedfile = pybedtools.BedTool(tempfile) # loading from tempfile
	bedfile = bedfile.sequence(fi = genome) # retrieving sequence within bedtools object
	fileout = list(SeqIO.parse(bedfile.seqfn, "fasta")) # creating a fasta object
	return fileout, bedfile.seqfn # 2 output type possible (SeqIO fasta object or Bedtool fasta file)
################################################
def reverse_complement(sequence):
	transtab = maketrans("ATGC","TACG")
	return sequence.upper().translate(transtab)[::-1]
################################################
def compute_TFFMscore(tffm,fasta, scoreTffmfile):
	sys.path.append("/home/304.3-STRUCTPDEV/lib/TFFM")
	import tffm_module
	from constants import TFFM_KIND
	hit_list_position=[]
	hit_list_score=[]
	tffm_first_order = tffm_module.tffm_from_xml(tffm, TFFM_KIND.FIRST_ORDER)
	with open(scoreTffmfile, "w") as OUT:
		for hit in tffm_first_order.scan_sequences(fasta, only_best=False):
			if str(hit)!="None":
				OUT.write(str(hit)+"\n")
################################################

def se_sp_matrix(posSet, negSet, pos_temp_filename, neg_temp_filename, genome, matrix, threhold):
	(waste, fastapos) = make_fasta(posSet, pos_temp_filename, genome) # getting fasta object
	(waste, fastaneg) = make_fasta(negSet, neg_temp_filename, genome) # getting fasta object
	if ".pfm" in matrix:
		import scores as sc # librairy created by Laura Gregoire
		sc.main(fastapos, matrix,".",False, "","") # calling scoring functions 
		sc.main(fastaneg, matrix, ".",False, "","") # calling scoring functions 
		list_pos = np.genfromtxt(fastapos.split("/")[-1] + ".scores", dtype = [('position','S75'), ('start',int), ('strand',float), ('score',float)]) # loading scoring results
		list_neg = np.genfromtxt(fastaneg.split("/")[-1] + ".scores", dtype = [('position','S75'), ('start',int), ('strand',float), ('score',float)]) # loading scoring results
	else:
		list_pos=compute_TFFMscore(matrix,fastapos, fastapos.split("/")[-1] + ".scores")
		list_neg=compute_TFFMscore(matrix,fastaneg, fastaneg.split("/")[-1] + ".scores")
		list_pos = np.genfromtxt(fastapos.split("/")[-1] + ".scores", dtype = [('position','S75'), ('start',int), ('end',int), ('strand',float), ('seq', 'S30'), ('mat','S4'), ('length',int),('score',float)]) # loading scoring results
		list_neg = np.genfromtxt(fastaneg.split("/")[-1] + ".scores", dtype = [('position','S75'), ('start',int), ('end',int), ('strand',float), ('seq', 'S30'), ('mat','S4'), ('length',int),('score',float)]) # loading scoring results
	se_list = []
	sp_list = []
	for thr in threhold: # for each threholds selected by user
		list_pos2 = np.unique([i for v,i in zip(list_pos['score'], list_pos['position']) if v >= float(thr)]) # unique position where score is superior to threhold
		length_pos = np.unique([i for i in list_pos['position']]) # unique positions
		list_neg2 = np.unique([i for v,i in zip(list_neg['score'], list_neg['position']) if v >= float(thr)]) # unique position where score is superior to threhold
		length_neg = np.unique([i for i in list_neg['position']]) # unique positions
		TruePos = list_pos2.size # number of sequences with score >= threhold and should
		FalseNeg = length_pos.size - list_pos2.size # number of sequences with score >= threhold and shouldn't
		FalsePos = list_neg2.size # number of sequences with score < threhold and shouldn't
		TrueNeg = length_neg.size - list_neg2.size # number of sequences with score > threhold and should
		se = float(TruePos) / float(TruePos + FalseNeg) # sensitivity
		spe = float(TrueNeg) / float(TrueNeg + FalsePos) # specificity
		se_list.append(se) # saving sensitivy for threhold considered 
		sp_list.append(1-spe) # saving specificity for threhold considered
	#positiveSet = np.flipud(np.sort(list_pos, order = 'score')) # sort by decreasing score
	#Set = np.unique(positiveSet['position'], return_index = True) # index of unique position
	#positiveSet = positiveSet[Set[1]]['score'] # best score of each position
	#negativeSet = np.flipud(np.sort(list_neg, order = 'score')) # sort by decreasing score
	#Set2 = np.unique(negativeSet['position'], return_index=True) # index of unique position
	#negativeSet = negativeSet[Set2[1]]['score'] # best score of each position
	os.remove(fastapos.split("/")[-1] + ".scores") # removing temp file 
	os.remove(fastaneg.split("/")[-1] + ".scores") # removing temp file
	return se_list, sp_list
################################################			   MAIN				################################################

results = {} # dictionnary for results to plot
limit_se = ""
limit_sp = ""
colors = ["b*","gv", "ms","c|","yx","kh"]
idcolor = 0
for elt in matrix:
	# calculating and saving sensitivity for each matrix
	(se_list, sp_list) = se_sp_matrix(posSet, negSet, pos_temp_filename, neg_temp_filename, genome, elt, threhold)
	results[str(elt)] = {}
	results[str(elt)]['se'] = se_list
	results[str(elt)]['sp'] = sp_list
	#results[str(elt)]['pos'] = positiveSet.tolist()
	#results[str(elt)]['neg'] = negativeSet.tolist()
	results[elt]['threhold'] = threhold
	results[elt]['color'] = "ro"

	#groups = ([1]*len(positiveSet.tolist())) + ([0]*len(negativeSet.tolist())) # attribuing groups 1 to posisitve set and 0 to negative set
	#scores = np.concatenate((positiveSet.tolist(), negativeSet.tolist())) # concatenate scores
	#fpr, tpr, _= sklearn.metrics.roc_curve(groups, scores, pos_label = 1) # computing ROC
	#results[str(elt)]['fpr'] = fpr # saving for later use
	#results[str(elt)]['tpr'] = tpr # saving for later use


print("="*20)
for elt in results.keys():
	print("results for "+elt+":")
	print("sensibilité :"+str(results[elt]['se']))
	print("spécificité :"+str(results[elt]['sp']))
	print("="*20)
try:
	os.remove(pos_temp_filename)
	os.remove(neg_temp_filename)
	#os.remove(pos_scoreTffmfile)
	#os.remove(neg_scoreTffmfile)
except:
	print("ok")


#python bin/thresholds_computation.py -g /home/304.6-RDF/data/tair10.fas -m /home/304.6-RDF/data/ARF2.pfm -d /home/304.6-RDF/Papier_ARF/Arf_Paper_fig1_maker/results/2018_02_12/ARF2_MP/Find_motifs/ARF2/testing_set/testing_set_pos.bed /home/304.6-RDF/Papier_ARF/Arf_Paper_fig1_maker/results/2018_02_12/ARF2_MP/Find_motifs/ARF2/testing_set/testing_set_1_neg.bed -t -8 -9 -10 -11 -12


#python bin/thresholds_computation.py -g /home/304.6-RDF/data/tair10.fas -t 0.4 0.35 0.3 -d /home/304.3-STRUCTPDEV/MADS/mapping_DAP/pipepline_rep/results/reanalysis/spacing/SEP3AG/sets/interdist_peaks_pos.bed /home/304.3-STRUCTPDEV/MADS/mapping_DAP/pipepline_rep/results/reanalysis/spacing/SEP3AG/sets/interdist_peaks_1_neg.bed -m /home/304.3-STRUCTPDEV/MADS/mapping_DAP/pipepline_rep/results/reanalysis/motifs/SEP3AG/SEP3AG_tffm.xml
