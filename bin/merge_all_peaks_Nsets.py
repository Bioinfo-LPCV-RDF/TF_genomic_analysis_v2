# -*- coding: utf-8 -*-
"""
merge_all_peaks_Nsets
=====================

DESCRIPTION
This script processes genomic intervals from a BED file to group them into contiguous regions (contigs) 
and merge overlapping peaks based on a specified overlap threshold. The merged peaks are written to an 
output file in BED format.
The script includes the following functionalities:
1. Parsing a BED file to group genomic intervals into contigs based on overlap criteria.
2. Merging contigs into peaks using simple or complex merging logic.
3. Writing the merged peaks to an output file.

USAGE
	python merge_all_peaks_Nsets.py -b <input_bedfile> -o <output_file> -po <overlap_percentage>

ARGUMENTS
	-b, --bedfile: Path to the sorted BED file containing genomic intervals. The file should be in 
				   BED format with columns: chromosome, start, stop, name, and optionally a score.
	-o, --output: Path to the output file where merged peaks will be written.
	-po, --percentageoverlap: Overlap percentage (default: 0.8) required for peaks to be considered 
							   as merged. The value should be between 0 and 1.

FUNCTIONS
	- define_contig(bedfile, overlap): Groups genomic intervals into contigs based on overlap criteria.
	- simple_case(a, b, overlap): Determines if two intervals can be merged using simple logic.
	- all_took(liste): Checks if all elements in a list are processed.
	- get_start(peak): Returns the start position of a peak.
	- complex_case(contigs, overlap): Handles complex merging scenarios for contigs.
	- define_merged(contigs, overlap, output): Merges contigs into peaks and writes them to an output file.
	- main(bedfile, overlap, output): Main function to process the BED file and generate merged peaks.
	- The script assumes the input BED file is sorted by chromosome and start position.
	- Overlap is calculated as the fraction of the smaller interval's length that overlaps with the other interval.
	- The script categorizes contigs into unique, simple, and complex types for merging.
	- The output file contains merged peaks in tab-delimited format: chromosome, start, end, name.

""" 
#!usr/bin/python
################################################                General intels                ################################################                
#description     : 
#author          : Jérémy Lucas
#date            : 18/01/2020
#version         :0.1
#usage           : [script name] -h
#python_version  :2.7.9
################################################                import packages                ################################################    
# library
#import os
#import metaseq
#import multiprocessing
#processes = multiprocessing.cpu_count()
#from pybedtools import BedTool
#import numpy as np
#from matplotlib import pyplot as plt
#from math import log10
#import matplotlib.cm as cm
#import matplotlib.colors as colors
#from matplotlib import gridspec

import argparse



################################################                Functions                ################################################
def define_contig(bedfile,overlap):
	"""
	Parses a BED file and groups genomic intervals into contiguous regions (contigs) 
	based on a specified overlap threshold.
	Args:
		bedfile (str): Path to the BED file containing genomic intervals. Each line 
					   in the file should represent an interval with at least 4 tab-separated 
					   columns: chromosome, start, end, and name. Additional columns are optional.
		overlap (float): Minimum overlap ratio (0 to 1) required to group two intervals 
						 into the same contig. The overlap is calculated as the fraction 
						 of the smaller interval's length that overlaps with the other interval.
	Returns:
		list: A list of contigs, where each contig is a list of intervals. Each interval 
			  is represented as a list containing its chromosome, start, end, name, 
			  and additional attributes if present in the input file.
	Notes:
		- Intervals are grouped into contigs if they are on the same chromosome and 
		  have sufficient overlap based on the specified threshold.
		- The function assumes the input BED file is sorted by chromosome and start position.
		- Overlap is calculated as the fraction of the smaller interval's length that 
		  overlaps with the other interval.
		"""
	last=[]
	contigs=[]
	potential=[]
	with open(bedfile) as IN:
		for line in IN:
			columns=line.split("\n")[0].split("\t")
			if len(columns)==4:
				current=[columns[0],int(columns[1]),int(columns[2]),columns[3],int((int(columns[2])-int(columns[1]))/2.0)]
			elif len(columns)==5:
				current=[columns[0],int(columns[1]),int(columns[2]),columns[3],int(columns[4])]
			else:
				current=[columns[0],int(columns[1]),int(columns[2]),columns[3],int(columns[4])]
				for elt in columns[5:]:
					current.append(elt)
			finish_contig=True
			if len(last)==0: # first line of the file
				last=current
				potential=[current]
				continue
			
			if current[0]==last[0]: # same chromosome
				if current[2]<=last[2]: # stop "current" before stop "last" => "current" contained in "last"
					potential.append(current)
					finish_contig=False # keep searching for peaks with overlap
				elif current[1]<=last[2]: # start "current" before stop "last" => partial overlap at least
					if (current[2]-current[1])>(last[2]-last[1]): # len(current) > len(last)
						if ((last[2]-current[1])/float(last[2]-last[1]))>=overlap: # overlap at least N% of the smallest peak
							potential.append(current)
							finish_contig=False # keep searching for peaks with overlap
					else: # len(current) <= len(last)
						if ((last[2]-current[1])/float(current[2]-current[1]))>=overlap: # overlap at least N% of the smallest peak
							potential.append(current)
							finish_contig=False # keep searching for peaks with overlap
			if finish_contig: # no more peaks with enough overlap
				contigs.append(potential)
				potential=[current] # last peak not in previous contig, create a new contig
			last=current # saving current peak as last before taking a new one
	contigs.append(potential) # saving last contig
	return contigs

def simple_case(a,b,overlap):
	distmax=(b[1]+b[4])-(a[1]+a[4])
	lenmin=min([b[2]-b[1],a[2]-a[1]])
	if distmax <= (((1-overlap)*2)*lenmin):
		return True,distmax
	else:
		return False,distmax

def all_took(liste):
	for elt in liste:
		if elt:
			continue
		else:
			return False
	return True

def get_start(peak):
	return peak[1]

def complex_case(contigs,overlap):
	new_peaks=[]
	took=[]
	for peak in contigs:
		took.append(False)
	continu=True
	i=-1
	while continu:
		print(str(contigs[i]))
		first=contigs[i][1]+contigs[i][4]
		lenmin=contigs[i][2]-contigs[i][1]
		newpeaks=first
		nb=1
		name=[contigs[i][3]]
		j=0
		for peaks in contigs:
			maximum=peaks[1]+peaks[4]
			lenpeak=peaks[2]-peaks[1]
			lenmin=min([lenmin,lenpeak])
			if abs(maximum-first)<=(((1-overlap)*2)*lenmin):
				newpeaks+=maximum
				nb+=1
				if peaks[3] not in name:
					name.append(peaks[3])
				took[j]=True
			j+=1
		name.sort()
		name = list(set(name))
		new_peaks.append([contigs[i][0],(newpeaks/nb)-(((1-(overlap))*2)*lenmin),(newpeaks/nb)+(((1-(overlap-0.1))*2)*lenmin),"&".join(name),newpeaks/nb])
		
		if i<0:
			i=0
		else:
			i+=1
		if all_took(took):
			continu=False
		elif i>0:
			i=took.index(False)
			print("results not 100% accurate")
			continu=False
	new_peaks.sort(key=get_start)
	return new_peaks

def define_merged(contigs,overlap,output):
	"""
	Processes genomic contigs to define merged peaks based on overlap criteria 
	and writes the resulting peaks to an output file.

	Args:
		contigs (list): A list of contigs, where each contig is a list of genomic 
						intervals. Each interval is represented as a list or tuple 
						with the following structure:
						[chromosome, start, end, name, score].
		overlap (float): A threshold value (0 to 1) that determines the overlap 
						 criteria for merging intervals.
		output (str): The file path where the merged peaks will be written.

	Returns:
		None

	The function categorizes contigs into three types:
		- Unique: Contigs with only one interval.
		- Simple: Contigs that can be merged based on the overlap criteria.
		- Complex: Contigs that require more advanced merging logic.

	The merged peaks are written to the output file in the following tab-delimited format:
		chromosome    start    end    name

	Additionally, the function prints a summary of the number of unique, simple, 
	and complex contigs processed, as well as the total number of contigs.

	Notes:
		- The function relies on two helper functions, `simple_case` and `complex_case`, 
		  which are not defined in the provided code.
		- The `simple_case` function is expected to handle simple merging logic, 
		  while `complex_case` handles more advanced merging scenarios.
	"""
	peaks=[]
	unique=0
	simple=0
	complexe=0
	for contig in contigs:
		if len(contig)==1:
			peaks.append(contig[0])
			unique+=1
		else:
			(result,dist)=simple_case(contig[0],contig[-1],overlap)
			if result:
				mean=0
				lenmin=contig[0][2]-contig[0][1]
				name=[]
				for elt in contig:
					mean+=elt[1]+elt[4]
					name.append(elt[3])
					if lenmin>(elt[2]-elt[1]):
						lenmin=(elt[2]-elt[1])
				name.sort()
				name = list(set(name))
				peaks.append([contig[0][0],(mean/len(contig))-(((1-(overlap))*2)*lenmin),(mean/len(contig))+(((1-(overlap-0.1))*2)*lenmin),"&".join(name),(mean/len(contig))])
				simple+=1
			else:
				print("=======")
				print(contig)
				print(dist)
				peaks+=complex_case(contig,overlap)
				complexe+=1
	with open(output,"w") as OUT:
		for elt in peaks:
			OUT.write(elt[0]+"\t"+str(int(elt[1]))+"\t"+str(int(elt[2]))+"\t"+elt[3]+"\n")
	print("=======")
	print("unique  : "+str(unique))
	print("multiple: "+str(complexe+simple))
	print("   simple  : "+str(simple))
	print("   complexe: "+str(complexe))
	print("total   : "+str(complexe+simple+unique))
	
################################################                MAIN                ################################################
def main(bedfile,overlap,output):
	contigs=define_contig(bedfile,overlap)
	#i=0
	#for elt in contigs:
		#i+=1
		#print(i)
		#for el2 in elt:
			#print("\t"+str(el2))
	define_merged(contigs,overlap,output)

################################################                Arguments                ################################################ 
if __name__ == "__main__": 
	parser = argparse.ArgumentParser() 
	
	parser.add_argument("--bedfile","-b",required=True, help="SORTED bedfile containing sets to compare. File should be in bed format (Chromosome, start, stop , name, maximum of heigt; name should be different between sets)")
	parser.add_argument("--output","-o",required=True)
	parser.add_argument("--percentageoverlap","-po",default =0.8,type=float,help="overlap percentage needed for peaks to be considered as merged")
	args = parser.parse_args()
	
	
	main(
		args.bedfile,
		args.percentageoverlap,
		args.output
		)
	
	#parser.add_argument("--bamdir","-b",required =True)
	#parser.add_argument("--peakfile","-p",required =True)
	#parser.add_argument("--resultdir","-r",required =True)
	#parser.add_argument("--centered","-c",required =True)
	#parser.add_argument("--windowssize","-w", type=int,required =True)
	#parser.add_argument("--orderwindow","-o",required =True)
##	parser.add_argument("--CFR","-c",default =[1,2],nargs="+",type=int)
	#parser.add_argument("--sort","-s", type=int, choices=[1,2,3,4],help = "1: sort on TF1; 2: sort on TF2; 3: sort on CFR",required =True)
	#parser.add_argument("--name","-n", nargs=2, help="name of the samples provided, in same order as in the peak file",required =True)
##	parser.add_argument("--all","-a",action='store_true', default= False)
	#
	#data_dir = args.bamdir
	#peak_file = args.peakfile
	#results = args.resultdir
	#sortBy = args.sort
	#list_name = args.name
	#windowssize = args.windowssize
	#order_window = args.orderwindow.split(":")
	#if args.centered == "True":
		#centered = True
	#else:
		#centered = False
	#main(data_dir, peak_file, results, sortBy, list_name, windowssize, centered,order_window)