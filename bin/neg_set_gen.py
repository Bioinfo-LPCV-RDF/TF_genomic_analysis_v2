#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Negative Set Generator for Genomic Analysis
============================================

DESCRIPTION:
This script generates negative sets of genomic regions for use in comparative
genomic analyses. It creates negative regions that match the characteristics of
positive regions (from ChIP-seq or similar experiments) in terms of GC content,
region type (exon, intron, promoter, etc.), and length distribution.

RELATED TO compute_NS

AUTHORS:
- Jérémy LUCAS (inspired by Arnaud STIGLIANI's work)
- Date: 16-01-2018

USAGE:
```bash
python neg_set_gen.py -pos <positive_regions> -fas <genome_fasta> -bed <genome_annotation>
                     [-of <output_prefix>] [-od <output_dir>] [-n <num_sets>]
                     [-r <random_seed>] [-GC <delta_GC>] [-l <limit_types>]
                     [-bs <bin_size>] [-st]
```

ARGUMENTS: 
-pos, --positive,    : BED file of positive regions (3 columns: chromosome, start, stop). Required.
-fas, --fasta        : Genome FASTA file.
-bed, --bedfile      : Genome annotation BED file (4 columns: chrom, start, stop, region_type). 
-of, --output_file   : Prefix for output files. Required.
-od, --output_dir    : Output directory.
-n, --amount         : Number of negative sets to generate. Default: 1.
-r, --rand           : Random seed. Default: 1.
-GC, --deltaGC       : Maximum allowed divergence of GC content between negative and positive sequences. Default: 0.03.
-l, --limittype      : Minimum number of sequences necessary to keep a region type. Default: 1000.
-bs, --bin_size      : Smallest bin size for genome partitioning. Default: 250.
-st, --simpletype    : Limit region types to simple origins (exon, intron, etc.). Flag.

EXAMPLE: 
#python /home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/bin/neg_set_gen.py -pos /home/312.6-Flo_Re/312.6.1-Commun/TF_benchmark/results/HB/ATHB5/Peakcalling/ChIPABA/ChIPABA_narrow10.bed -of ChIPABA -od /nobackup//NSG/ChIPABA -fas /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas -bed /home/312.6-Flo_Re/312.6.1-Commun/data/tair10_prepJL.bed -n 1 -r 1254 -GC 0.03 -l 6000 -bs 250

OUTPUTS: 
- <output_prefix>_pos.bed       : Processed positive regions
- <output_prefix>_<n>_neg.bed   : Negative regions for each requested set (n=1,2,...)

NOTES: 
- Region types can be: intergenic, promoter, 5UTR, 3UTR, intron, exon
- The bin_size parameter should be approximately the mean length of positive regions + 50
- For small genomes or region sets, the limittype parameter may need adjustment
"""

################################################			 import packages				################################################			 
import argparse # version 1.1
from Bio import SeqIO # version 1.70
import pybedtools # version 0.7.10
import random 
import numpy as np # version 1.8.2
import time 
import os
import sys
import tempfile
start = time.time()


################################################			 Parsing arguments				################################################
parser = argparse.ArgumentParser(description = "Create appropriate negative set using data from Genome and positive set (type of region, GC_content, length)") 
########### 
parser.add_argument("--positive","-pos",help = "Bed file of your positives regions, it may contain 3 columns tab separated (chromosome, start, stop)", required = True)
parser.add_argument("--fasta","-fas",default = "../data/tair10.fas", help = "fasta of your genome, header for each chromosome should match those specified in your positive regions", required = True)
parser.add_argument("--bedfile","-bed",default = "../data/Genome_info.bed", help = "Bed file of your genome, it may contain 4 columns tab separated (chromosome, start, stop, type_of_regions). types of regions accepted are : \"intergenic, promoter, 5UTR, 3UTR, intron, exon\". All non accepted regions will be considered as \"intergenic\"", required = True)
parser.add_argument("--bin_size","-bs",type = int,default = 250, help = "smallest bin size for your genome. [default: 250], a good choice may be the mean length of your positive regions + 50. If your choice is too small for some or all regions, this bin will be multiplied accordingly")
parser.add_argument("--output_file","-of",help = "prefix of your output files (negative and positive set files)", required = True)
parser.add_argument("--amount","-n",type = int,default = 1, help = "number of negative sets to generate [default: 1]")
parser.add_argument("--output_dir","-od",default = "./", help = "result directory. [default: ./]")
parser.add_argument("--rand","-r",type = int,default = 1, help = "random seed [default: 1]")
parser.add_argument("--deltaGC","-GC",type = float,default = 0.03, help = "maximum authorized divergence of GC content betwen negative and positive sequences")
parser.add_argument("--limittype","-l",type = int,default = 1000, help = "minimum number of sequences necessary to keep region type.")
parser.add_argument("--simpletype","-st",action="store_true",default = False, help = "limit the region type to simple origins (Exon,intron, etc)")
########### 
args = parser.parse_args()
amount_of_negative_set = args.amount
name_neg_file = args.output_file
name_pos_file = args.positive
output_dir = args.output_dir
genome_bin_size = args.bin_size
name_bed_file = args.bedfile
name_fasta_file = args.fasta
seed_rand = args.rand
limit_types = args.limittype
deltaGC = args.deltaGC
simpletype=args.simpletype
genome_regions_library = {}
########### 
if not output_dir.endswith("/"):
	output_dir += "/"
if simpletype:
	print("simpletype is active, limit type is ignored")
################################################			 creating random tempfiles				################################################
temp_var = tempfile.NamedTemporaryFile(dir = output_dir)
positive_temp_filename = temp_var.name
temp_var.close()
temp_var = tempfile.NamedTemporaryFile(dir = output_dir)
positive_extended_temp_filename = temp_var.name
temp_var.close()
temp_var = tempfile.NamedTemporaryFile(dir = output_dir)
general_temp_filename = temp_var.name
temp_var.close()
################################################			 Functions				################################################
'''
allow the truncature of each feature taken from an intersect of bedtools objects
from a intersect, create a bedtools object of 4 columns (chromosome, start, stop, annotation)
'''
def trunc_feature(feature):
	new_feature = feature[0:3] # retrieve chromosome, start, stop
	new_feature.append(feature[6]) # retrieve annotation
	new_feature.append(feature[7]) # retrieve overlap with annotation
	new_feature = pybedtools.create_interval_from_list(new_feature) # recreate interval to keep integrity of bedtools object
	return new_feature
################################################
'''
compute GC content of a DNA region and filter out sequences which contain at least one "N" (or anything not in the alphabet : ATGC)
'''
def GC_content(chrom,start,stop,precision):
	global Fasta_tair
	try:
		local_fasta = Fasta_tair[chrom].seq[start-1:stop].upper() # retrieving DNA seq for 1 region
		if len(local_fasta) != (stop-start)+1: # region not entirely in chromosome
			return "" # stops research
		# computing GC content
		G = local_fasta.count('G')
		C = local_fasta.count('C')
		GC = float(G+C)/(len(local_fasta))
		# trying to figure out if ATGC are the only letters in the sequence alphabet
		N = (len(local_fasta)) -(local_fasta.count('A')+local_fasta.count('T')+local_fasta.count('G')+local_fasta.count('C'))
		if N == 0:
			return round(GC,precision)
		else: # if there is more letters in the alphabet than ATGC
			return "" # regions is discarded
	except:
		print("\tunable to get region:")
		print("\t"+str(chrom)+"\t"+str(start)+"\t"+str(stop))
		print(str(local_fasta))
		print(str(G),str(C),str(GC))
		return "" # stops research
################################################
"""
save values with same key in a string, with a ";" as separator
"""
def create_dict(dicts,key,value):
	if key not in dicts.keys():
		dicts[key]=str(value)
	else:
		dicts[key]+=";"+str(value)
	return dicts
################################################
'''
merge positive regions, created after intersected with genome information bedtools object, which share same location and create a suitable bedtools object for following operations
'''
def merge_features(filename):
	positive_dataset = np.genfromtxt(filename,dtype = [('chrom','U20'),('start',int),('stop',int),('type','U300')]) # loading positive set
	appendice = []
	keys = set()
	dicts = {}
	# create a new numpy array with no feature on same position
	appendice = [keys.add(str(sublist["chrom"])+"-"+str(sublist["start"])+"-"+str(sublist["stop"])) or sublist 
				  for sublist in positive_dataset 
					  if str(sublist["chrom"])+"-"+str(sublist["start"])+"-"+str(sublist["stop"]) not in keys]
	del keys 
	trash =[]
	# create a dictionnary with every annotation found for each position
	trash = [create_dict(dicts,str(sublist["chrom"])+"-"+str(sublist["start"])+"-"+str(sublist["stop"]),sublist["type"])
			  for sublist in positive_dataset]
	# merge unique position numpy array and annotations
	appendice= [[str(sublist[0]),int(sublist[1]),int(sublist[2]),dicts[str(sublist[0])+"-"+str(sublist[1])+"-"+str(sublist[2])]] or sublist
						 for sublist in appendice 
						 if str(sublist[0])+"-"+str(sublist[1])+"-"+str(sublist[2]) in dicts.keys()]
	new_positive_dataset = pybedtools.BedTool(appendice) # creating a bedtools object
	new_positive_dataset = new_positive_dataset.saveas(filename+"_2") # saving Object
#	test = new_positive_dataset.intersect(new_positive_dataset, wo = True, f = 1.00)
################################################
'''
allow the truncature of each feature taken from an intersect of bedtools objects
filtering annotation to keep only combination of types allowed [intergenic, promoter, 5UTR, 3UTR, intron, exon]
all unauthorized combination are discarded and, if no annotation remains, replace by default by intergenic
store GC content and length of each feature
'''
def trunc_feature_filter(feature):
	new_list = []
	# detection of multiple annotation for the feature
	word_list_temp = feature.name.split(";")
	word_list=[]
	for elt in word_list_temp:
		word_list.append(elt.split("/")[0])
	# filtering out RNA annotation and pseudogenic_transcript, keeping only types from this list
	word_list = list(filter(lambda x: x in ["exon","intron","intergenic","5UTR","3UTR","promoter"],word_list))
	for word in word_list: # truncates annotation with name of gene attached after a "/"
		if word not in new_list: # avoid duplicates in annotation
			new_list.append(word)
	if len(new_list) == 0: # no annotation left after filtering
		new_list = ["intergenic"] # default annotation
	new_list.sort() # sorting annotation to avoid multiple times same types in different orders
	feature.name = ";".join(new_list)
	GC = str(GC_content(feature.chrom,feature.start,feature.stop,1)) # call function for GC content storage
	length = feature.stop-feature.start # length of feature storage
	if GC  != "" and GC  != 2: # if GC content calculated (feature not out of chromosomes)
		new_feature = feature
		new_feature.append(GC)
		new_feature.append(length)
		return new_feature
################################################
'''
Regroups steps for treatement of bedtools object
'''
def prepare_bedtools_object(temp_save,simpletype):
	global info_regions
	Object = pybedtools.BedTool(temp_save) # loading bedtools object from tempfile
	if simpletype:
		Object = Object.intersect(info_regions, wo = True, f = 0.10,nonamecheck=True).saveas(temp_save) # intersect with annotated bedtools object to transfert annotations
		test_sortinter = np.genfromtxt(temp_save,dtype = [('chrom','U20'),('start',int),('stop',int),('chrom2','U20'),('start2',int),('stop2',int),('type','U300'),('intersect',int)])
		a = np.sort(test_sortinter, order=['chrom','start','intersect'])
		result = a[np.concatenate((a['start'][:-1] != a['start'][1:], np.array([True])))]
		liste = result.tolist()
		Object = pybedtools.BedTool(liste)
		
	else:
		Object = Object.intersect(info_regions, wo = True, f = 0.10,nonamecheck=True)
	######
	# Object = Object.saveas(temp_save+"bis")
	# Object = pybedtools.BedTool(temp_save+"bis")
	# sys.exit()
	Object = Object.each(trunc_feature).saveas(temp_save) # reformat bedtools object, sorting and saving it to tempfile
	# print(Object) # new length of positive region, after intersect with genome info
	merge_features(temp_save)
	
	Object = pybedtools.BedTool(temp_save+"_2") # loading bedtools object from tempfile	
	# sys.exit()
	
	Object = Object.each(trunc_feature_filter) # integrating length, GC content and filtering annotations
	Object = Object.saveas(temp_save)
	Object = pybedtools.BedTool(temp_save)
################################################
'''
OUTPUT function for positives regions
'''
def print_region(positions_pos_regions):
	global output_dir
	global name_neg_file
	with open(output_dir+name_neg_file+"_pos.bed", "w") as OUT:
#		 map(lambda x: OUT.write("\t".join(x[0:3])+"\n"),positions_pos_regions)
		for feature in positions_pos_regions:
			line = "\t".join(feature[0:3])
			OUT.write(line+"\n")
	OUT.close()
################################################
'''
extend length of feature by 2 (1bp on each side) and return the figure
'''
def extend_pos_reg(feature):
	feature.start -= 1
	feature.stop += 1
	return feature
################################################
'''
Regroups functions for the treatment of positives regions
'''
def treat_positive_region(name_pos_file,name_fasta_file,name_bed_file,simpletype):
	global info_regions
	global positive_temp_filename
	global positive_extended_temp_filename
	print("- Importing positive regions")
	positions_pos_regions = np.genfromtxt(name_pos_file,dtype = [('chrom','U20'),('start',int),('stop',int)])
	np.savetxt(positive_temp_filename,X=positions_pos_regions,delimiter="\t",fmt="%s") # save bed file in output directory
	positions_pos_regions = pybedtools.BedTool(positive_temp_filename)# loading bedtools object from input file
	print("- Annotating positives regions")
	prepare_bedtools_object(positive_temp_filename,simpletype) # working on bedtools object to annotate region, add GC content and length
	positions_pos_regions = pybedtools.BedTool(positive_temp_filename)# loading bedtools object from tempfile
	print("\tafter filtering invalid sequences "+str(len(positions_pos_regions))) # new length of positive region, without region containing Ns
	print_region(positions_pos_regions) # saving new positive set to output file
	Nmin = min(map(lambda x: int(x[2])-int(x[1]),positions_pos_regions)) # retrieving minimum length of positive region
	Nmax = max(map(lambda x: int(x[2])-int(x[1]),positions_pos_regions)) # retrieving maximum length of positive region
	temp_pos_region = positions_pos_regions.each(extend_pos_reg).saveas(positive_extended_temp_filename) # create a bed file with extended feature positive set. This will help to remove those from genome library
	return positive_temp_filename, Nmin, Nmax
################################################
'''
allow the truncature of each feature taken from an intersect of bedtools objects
from a intersect, create a bedtools object of 4 columns (chromosome, start, stop, annotation)
'''
def trunc_feature_genome(feature):
	new_feature = feature[0:3] # retrieve chromosome, start, stop
	new_feature.append(feature[7]) # retrieve annotation
	new_feature = pybedtools.create_interval_from_list(new_feature) # recreate interval to keep integrity of bedtools object
	return new_feature
################################################
'''
Regroups steps for treatement of bedtools object containing whole genome data
'''
def prepare_bedtools_object_genome(temp_save,simpletype):
	global info_regions
	Object = pybedtools.BedTool(temp_save) # loading bedtools object from tempfile
	# print(Object.head(15))
	#Object = Object.intersect(info_regions, wo = True,f = 0.10,nonamecheck=True)#.saveas(temp_save) # intersect with annotated bedtools object to transfert annotations
	if simpletype:
		Object = Object.intersect(info_regions, wo = True, f = 0.10,nonamecheck=True).saveas(temp_save) # intersect with annotated bedtools object to transfert annotations
		test_sortinter = np.genfromtxt(temp_save,dtype = [('chrom','S20'),('start',int),('stop',int),('length',int),('chrom2','S20'),('start2',int),('stop2',int),('type','S300'),('intersect',int)])
		a = np.sort(test_sortinter, order=['chrom','start','intersect'])
		result = a[np.concatenate((a['start'][:-1] != a['start'][1:], np.array([True])))]
		liste = result.tolist()
		Object = pybedtools.BedTool(liste)
	else:
		Object = Object.intersect(info_regions, wo = True, f = 0.10,nonamecheck=True)
	
	
	Object = Object.each(trunc_feature_genome).saveas(temp_save) # reformat bedtools object, sorting and saving it to tempfile
	Object = pybedtools.BedTool(temp_save) # loading bedtools object from tempfile
	
	length = Object[0].stop-Object[0].start # get a length
	Object = Object.merge(c = 4,d = -length, o = "collapse", delim = ";") # collapse feature using length
	Object = Object.each(trunc_feature_filter) # integrating length, GC content and filtering annotations
	Object.saveas(temp_save)
	Object = pybedtools.BedTool(temp_save)
################################################
'''
cutting genome in bins of different length
'''
def cut_this_genome_in_bin(Nmin, Nmax,name_fasta_file,simpletype):
	global Fasta_tair
	global genome_bin_size
	global genome_regions_library
	global general_temp_filename
	# calculating number of bins we need
	imax = (Nmax//genome_bin_size) # maximum of times bin_sizes fits in maximum length of positive regions
	imin = (Nmin//genome_bin_size) # maximum of times bin_sizes fits in minimum length of positive regions
	rmax = Nmax-(imax*genome_bin_size) # remains after max, to consider all positive regions
	rmin = Nmin-(imin*genome_bin_size) # remains after max, to consider all positive regions
	if rmax >= -25:
		imax += 1
	if imin == 0 or rmin >=-25: # in case Nmin is to small
		imin += 1
	h = imin
	genome_bin = {}
	bin_size_list=[]
	while(h <= imax): # for each bin size
			length_reg = h*genome_bin_size # length of current bin
			bin_size_list.append(length_reg)
			print("\tworking on bin "+str(length_reg))
			h += 1
			chrom = []
			Save_cutted_genome = "/".join(name_fasta_file.split("/")[:-1])+"/"+name_fasta_file.split("/")[-1].split(".")[0]+"_"+str(length_reg)+".bed" # filename for saved cutted genome libraries in fasta directory
			if not os.path.isfile(Save_cutted_genome): # checking if saved cutted genome library exist
				for elt in Fasta_tair.keys():#  for each chromosome
					chrom.append("".join(list(map(lambda x: Fasta_tair[elt].id+"\t"+str(x)+"\t"+str(x+length_reg)+"\t"+str(length_reg)+"\n",range(1, len(Fasta_tair[elt].seq),length_reg+1))))) # make a list of bins on the entire chromosome
				genome = pybedtools.BedTool("".join(chrom), from_string = True).saveas(general_temp_filename+"_"+str(length_reg)) # creating a bedtools object from list and saving it
				prepare_bedtools_object_genome(general_temp_filename+"_"+str(length_reg),simpletype)
				genome = pybedtools.BedTool(general_temp_filename+"_"+str(length_reg)) # loading bedtools object from tempfile
				genome_bin[length_reg] = general_temp_filename+"_"+str(length_reg)
				genome.saveas(Save_cutted_genome) # saving for later use
			else:
				print("\t loading library from saved file: "+Save_cutted_genome)
				genome_bin[length_reg] = Save_cutted_genome
	return genome_bin,bin_size_list
################################################
def lambda_filter(x,tokeep):
	r=""
	y = x.split(";")
	if len(y) == 0:
		r = "intergenic"
	elif len(y) == 1:
		if y not in tokeep:
			r = "intergenic"
	elif len(y) >= 2:
		for elt in ["intergenic","3UTR","5UTR","exon","intron","promoter"]:
			if elt in y:
				w = x.split(";")
				w.remove(elt)
				if ";".join(w) in tokeep:
					r = ";".join(w)
					break
		if r == "":
			for elt in ["intergenic","3UTR","5UTR","exon","intron"]:
				if elt in y:
					w = x.split(";")
					w.remove(elt)
					r = ";".join(w)
					break
	return r
################################################
"""
for each bin size, features overlapping with positive set (length extended by two) are removed.
Types for genome library and positive set are filtered to limit small sets (harder to find suitables sequences) and faster execution.
"""
def filter_region(genome_regions_library, list_pos, limit_types, positive_extended_temp_filename, general_temp_filename,simpletype):
	for bin_length in genome_regions_library.keys(): # iterating on bin size
		temp_genome = pybedtools.BedTool(genome_regions_library[bin_length]) # loading genome library as bedtool object
		temp_pos_region = pybedtools.BedTool(positive_extended_temp_filename) # loading extended positive set as bedtool object
		temp_genome = temp_genome.subtract(temp_pos_region, A = True,nonamecheck=True).saveas(general_temp_filename+"_"+str(bin_length)) # removing positive set features from genome library
		genome_regions_library[bin_length] = np.genfromtxt(general_temp_filename+"_"+str(bin_length),dtype = [('chrom','U25'),('start',int),('stop',int),('type','U300'),('GC',float),('length',int)]) # loading genome library
		list_type_reg = np.unique([v for i,v in enumerate(genome_regions_library[bin_length]['type'])]) # capturing names for all detected types
		tokeep = []
		for elt in list_type_reg: # iterating on each detected types (single or composite)
			origin_genome = [i for i,v in enumerate(genome_regions_library[bin_length]['type']) if elt == v] # capturing number of sequences corresponding to type
			if simpletype :
				tokeep.append(elt)
			else:
				if len(origin_genome) >= limit_types and len(elt.split(";")) <= 3: # check if count if over "limit" and remove too composite types (>3) 
					tokeep.append(elt)
		origin_genome = [i for i,v in enumerate(genome_regions_library[bin_length]['type']) if v in tokeep] # filtering out type which are not kept
		genome_regions_library[bin_length] = genome_regions_library[bin_length][origin_genome] # saving results
		stop = False
		i=0
		while not stop :
			"""
			while there is differences between types kept and types in the positives set, 1 types from composites sites are removed, UTR first and then intergenic.
			"""
			types_pos=np.unique([v for i,v in enumerate(list_pos['type'])]).tolist()
			if types_pos == tokeep :
				stop = True
			else:
				stop = True
				for elt in types_pos :
					if elt not in tokeep:   
						stop = False
						break
				if not stop:
					list_pos['type'] = [ lambda_filter(x,tokeep) for x in list_pos['type']]
			i+=1
			if i>=20:
				sys.exit() # safeguard for no infinite loop
	return genome_regions_library,list_pos
################################################
'''
delete bins used for negative set.
If the genome have been cut in different bin size, every bins of every size that overlap the corresponding selected bin will be deleted.

example with 4 differents bin size (250,500,750,1000)
considering those regions :
bins 250 :	250	 |	500	|	750	|	1000	|	1250	|
bins 500 :		500		|		1000		|		1500
bins 750 :		750		|			 1500		 
bins 1000:						1000					 |			 2000
if region 500 from bins 250 is selected as negative set, 
it's needed we delete overlapping bins of other sizes (region 500 from bins 500, region 750 from bins 750 etc)

to do so safely (aka not delete a region for nothing or keep a region overlapping an used region), we use this formula:
length_reg*(position of region ) // genome_bin_size *i = first region to delete

if length_reg*(position of region +1) % genome_bin_size *i > 0:
	(length_reg*(position of region +1) // genome_bin_size *i) = last region to delete
else:
	(length_reg*(position of region +1) // genome_bin_size *i)-1 = last region to delete
'''
def delete_from_library(positions_to_delete,imin,imax,genome_bin_size,length_reg):
	global temp_genome_lib
	to_del = {}
	i = imin
	while i <= imax: # preparing temp storage of all regions to delete (selected and overlapping)
		to_del[genome_bin_size*i] = []
		i += 1
	for position_to_delete in positions_to_delete: # for each region selected as negative set
		i = imin
		while i <= imax: # for each bin size created for the genome
			# applying formulas
			start_bin = (length_reg*position_to_delete)//(genome_bin_size*i) 
			if (length_reg*(position_to_delete+1))%(genome_bin_size*i) != 0:
				stop_bin = (length_reg*(position_to_delete+1))//(genome_bin_size*i)
			else:
				stop_bin = ((length_reg*(position_to_delete+1))//(genome_bin_size*i))-1
			j = start_bin
			while j <= stop_bin: # saving all regions to delete in temp storage
				to_del[genome_bin_size*i].append(int(j))
				j += 1
			i += 1
	for key,value in to_del.items(): # deleting
		# print(key)
		temp_genome_lib[key] = np.delete(temp_genome_lib[key],to_del[key])
################################################
'''
compute search of a negative sequence for each positive seqs given
search for identical types of region and same GC content (+-3%)
'''
def extract_neg_region(length_reg, origin_GC_genome, region_pos_considered, origin_GC_genome_approx, genome_bin_size, imin, imax, deltaGC, origin_GC_genome_closecall):
	global temp_genome_lib
	#origin_GC_genome_closecall = []
	list_neg = []
	positions_to_delete = []
	savegenome=origin_GC_genome_closecall
	initial_length = [len(region_pos_considered),len(origin_GC_genome),len(origin_GC_genome_closecall),len(origin_GC_genome_approx)] # save values before starting calculation to report errors
	#print("number of positive regions for this category"+"\t"+str(initial_length[0])+"\n"+"number of potential regions for this category"+"\t"+str(initial_length[1])+"\n"+"number of potential regions close to this category"+"\t"+str(initial_length[3])+" "+str(initial_length[2])+"\n")
	for region_pos in region_pos_considered:
		search_continue = True
		i = 0
		origin_GC_continue = True
		origin_GC_approx_continue = True
		origin_GC_genome_closecall_continue = True
		delta_i = deltaGC #max divergence of GC content between sequences
		tried_approx=[]
		tried_closecall=[]
		while search_continue: # while no sequence found
			 
			if origin_GC_genome.size != 0 and origin_GC_continue: # while there is non used sequence with same GC and type, use them
				#print("using orignal")
				random_choice = random.randint(0,len(origin_GC_genome)-1) # random choice of 1 seq in set
				position_region_neg = int(origin_GC_genome[random_choice])
				region_neg = temp_genome_lib[length_reg][position_region_neg] 
			elif len(origin_GC_genome_approx) != 0 and origin_GC_approx_continue : # regions with same type and mean GC close
				#print("Using genome_approx")
				random_choice = random.randint(0,len(origin_GC_genome_approx)-1)# random choice of 1 seq in set
				position_region_neg = int(origin_GC_genome_approx[random_choice])
				region_neg = temp_genome_lib[length_reg][position_region_neg]
				tried_approx.append(random_choice)
			elif len(origin_GC_genome_closecall) != 0 and origin_GC_genome_closecall_continue:
				# print("=====Using genome_closecall")
				random_choice = random.randint(0,len(origin_GC_genome_closecall)-1)# random choice of 1 seq in set
				position_region_neg = int(origin_GC_genome_closecall[random_choice])
				region_neg = temp_genome_lib[length_reg][position_region_neg]
				tried_closecall.append(int(origin_GC_genome_closecall[random_choice]))
			elif origin_GC_genome_closecall_continue == False:
				print("Warning - tried everything, not successfully, selecting random sequence")
				print(str(tried_closecall[0]))
				random_choice = random.randint(0,len(savegenome)-1)
				position_region_neg = int(savegenome[random_choice])
				region_neg = temp_genome_lib[length_reg][position_region_neg]
			else:
				print("ERROR") # no sequence left in library, should not be printed
				print("number of positive regions for this category"+"\t"+str(initial_length[0])+"\n"
					"number of potential regions for this category"+"\t"+str(initial_length[1])+"\n"
					"number of potential regions close to this category"+"\t"+str(initial_length[3])+" "+str(initial_length[2])+"\n")
				print("cannot find negative region for: "+"\t"+region_pos['chrom']+"\t"+str(region_pos['start'])+"\t"+str(region_pos['stop']))
				print("try with a deltaGC higher")
				print("tried "+str(len(tried_approx))+" "+str(len(tried_closecall)))
				delete_from_library(positions_to_delete,imin,imax,genome_bin_size,length_reg)
				return list_neg
			i += 1 
			delta_length = region_neg['length']-region_pos['length']
			GC_pos = GC_content(region_pos['chrom'],region_pos['start'],region_pos['stop'],2) # precise GC content for pos set
			h = 0
			delta = 0.01 # initial authorized divergence of GC
			while h <= delta_length:
				if origin_GC_genome_closecall_continue!=False:
					GC_neg = GC_content(region_neg['chrom'],region_neg['start']+h,region_neg['stop']-(delta_length-h),2)
				else:
					GC_neg=GC_pos
				h += 1
				if GC_neg != "" and GC_pos+delta >= GC_neg >= GC_pos-delta:
					#if ((len(origin_GC_genome) != 0 and origin_GC_continue) and origin_GC_genome[random_choice] not in positions_to_delete) or origin_GC_genome_approx[random_choice] not in positions_to_delete:
					list_neg.append(region_neg['chrom']+"\t"+str(region_neg['start']+h)+"\t"+str(region_neg['stop']-(delta_length-h)))
					if len(origin_GC_genome) != 0 and origin_GC_continue:
						if origin_GC_genome[random_choice] not in positions_to_delete:
							positions_to_delete.append(origin_GC_genome[random_choice])
							origin_GC_genome = np.delete(origin_GC_genome,random_choice)
					elif len(origin_GC_genome_approx) != 0 and origin_GC_approx_continue:
						if origin_GC_genome_approx[random_choice] not in positions_to_delete:
							positions_to_delete.append(origin_GC_genome_approx[random_choice])
							origin_GC_genome_approx = np.delete(origin_GC_genome_approx,random_choice)
					elif origin_GC_genome_closecall[random_choice] not in positions_to_delete:
						positions_to_delete.append(origin_GC_genome_closecall[random_choice])
						origin_GC_genome_closecall = np.delete(origin_GC_genome_closecall,random_choice)
					elif savegenome[random_choice] not in positions_to_delete:
						positions_to_delete.append(savegenome[random_choice])
						origin_GC_genome_closecall = np.delete(savegenome,random_choice)
					else:
						print(random_choice)
						print(origin_GC_genome_closecall[random_choice])
						if  origin_GC_genome_closecall[random_choice] in positions_to_delete:
							print("selected sequence already in positions to delete")
						print("initial: "+str(origin_GC_continue))
						print("approx: "+str(origin_GC_approx_continue))
						print("closecall: "+str(origin_GC_genome_closecall_continue))
						print("ERROR, sequence not found ERROR")
						print("number of positive regions for this category"+"\t"+str(initial_length[0])+"\n"+"number of potential regions for this category"+"\t"+str(initial_length[1])+"\n"+"number of potential regions close to this category"+"\t"+str(initial_length[3])+" "+str(initial_length[2])+"\n")
						quit()
					search_continue = False
					h = delta_length+1
				if h>delta_length and delta <= delta_i and not GC_pos+delta >= GC_neg >= GC_pos-delta:
					delta += 0.01
					h = 0
				if h > delta_length and delta > delta_i and not GC_pos+delta >= GC_neg >= GC_pos-delta and deltaGC*2 > delta_i and not origin_GC_continue and not origin_GC_approx_continue and len(origin_GC_genome_closecall) != 0 and len(origin_GC_genome_closecall) < i-1 and origin_GC_genome_closecall_continue :
				 # and (len(tried_approx) >= len(origin_GC_genome_approx) or not origin_GC_approx_continue) and (len(tried_closecall) >= len(origin_GC_genome_closecall) or len(origin_GC_genome_closecall)==0) and deltaGC*2 > delta_i:
					#print("tried_approx: "+str(len(tried_approx)))
					#print("tried_closecall: "+str(len(tried_closecall)))
					#print(search_continue)
					tried_approx=[]
					tried_closecall=[]
					origin_GC_continue = True
					origin_GC_approx_continue = True
					origin_GC_genome_closecall_continue = True
					print("WARNING: cannot find suitable negative sequence, changing deltaGC from "+str(delta_i)+" to "+str(round(delta_i+0.02,3)))
					delta_i = round(delta_i+0.02,3)
					i=0
				if h > delta_length and delta > delta_i and not GC_pos+delta >= GC_neg >= GC_pos-delta :
					if len(origin_GC_genome) < i-1 and origin_GC_continue:
						origin_GC_continue = False
						i=0
						#print("origin_GC ended")
					if len(origin_GC_genome_approx) != 0 and len(origin_GC_genome_approx) < i-1 and origin_GC_approx_continue:
						#print("end: "+str(len(origin_GC_genome_approx))+" "+str(len(tried_approx)))
						origin_GC_approx_continue = False
						i=0
						#print("origin_GC_approx ended")
					if len(origin_GC_genome_closecall) != 0 and len(origin_GC_genome_closecall) < i-1 and origin_GC_genome_closecall_continue:
						origin_GC_genome_closecall_continue = False
						i=0
						#print("origin_GC_closecall ended")
	delete_from_library(positions_to_delete,imin,imax,genome_bin_size,length_reg)
	return list_neg
################################################
'''

'''
def fit_negative_set(temp_genome_lib,list_pos, genome_bin_size,deltaGC):
#	list_pos = np.genfromtxt(positions_pos_regions_filename,dtype = [('chrom','S20'),('start',int),('stop',int),('type','S300'),('GC',float),('length',int)])
	list_neg = []
	# calculating number of bins used
	imax = (Nmax//genome_bin_size) # maximum of times bin_sizes fits in maximum length of positive regions
	imin = (Nmin//genome_bin_size) # maximum of times bin_sizes fits in minimum length of positive regions
	rmax = Nmax-(imax*genome_bin_size) # remains after max, to consider all positive regions
	rmin = Nmin-(imin*genome_bin_size) # remains after max, to consider all positive regions
	if rmax >= -25:
		imax += 1
	if imin == 0 or rmin >=-25: # in case Nmin is to small
		imin += 1
	leng = imin
	while(leng <= imax): # for each bin size
		length_reg = leng*genome_bin_size 
		list_type_reg = np.unique([v for i,v in enumerate(list_pos['type'])]).tolist() # considering all types in positives regions
		print("\tworking on bin "+str(length_reg))
		leng += 1
		h = 0.0
		length_pos = [i for i,v in enumerate(list_pos['length']) if length_reg == (((v//genome_bin_size)+1)*genome_bin_size)] # retrieving all positive region of this length
		while h <= 1: # for each mean of GC content
			GC_searched = round(h,1)
			min_gc_display = round((h-0.05)*100,0)
			max_gc_display = round((h+0.05)*100,0)
			if min_gc_display < 0.0:
				min_gc_display = 0.0
			if max_gc_display > 100:
				max_gc_display = 100.0 
			h += 0.1
			GC_pos = np.where(list_pos['GC'] == GC_searched)[0].tolist() # retrieving all positive regions of this mean GC content
			GC_pos = np.intersect1d(GC_pos,length_pos) # intersecting result to get regions of considerated length and type
			print("\t  working on GC ["+str(min_gc_display)+"%;"+str(max_gc_display)+"%[, "+str(len(GC_pos))+" sequences to process")
			for elt in list_type_reg: # for each type considered
				origin_GC_genome_approx = []
				origin_GC_genome_closecall = []
				origin_pos = [i for i,v in enumerate(list_pos['type']) if elt == v] # retrieving all positive regions of this type
				origin_GC_pos = np.intersect1d(GC_pos,origin_pos) # intersecting result to get regions of considerated length, mean GC content and type
				if len(origin_GC_pos)  != 0: # if there is positive regions to treat
					origin_genome = [i for i,v in enumerate(temp_genome_lib[length_reg]['type']) if elt == v] # retrieving all genomic regions of this type and length
					GC_genome = np.where(temp_genome_lib[length_reg]['GC'] == GC_searched)[0].tolist() # retrieving all genomic regions of this mean GC content and length
					origin_GC_genome = np.intersect1d(GC_genome,origin_genome) # intersecting result to get regions of considerated length, mean GC content and type
					if origin_GC_genome.size<origin_GC_pos.size+origin_GC_pos.size/2:
						i = 1
						GC_genome = []
						while (origin_GC_genome.size+len(origin_GC_genome_approx))<=round(origin_GC_pos.size+origin_GC_pos.size/2,0) and i <= 4:
							new_limit_GC = i*0.1
							i += 1
							new_limit_inf= max(round(GC_searched-new_limit_GC,1),0.0)
							new_limit_sup= min(round(GC_searched+new_limit_GC,1), 1.0)
							#print(new_limit_inf)
							#print(new_limit_sup)
							try:
								
								GC_genome = np.concatenate((GC_genome,np.where(temp_genome_lib[length_reg]['GC'] == new_limit_inf)[0].tolist()))
								GC_genome = np.concatenate((GC_genome,np.where(temp_genome_lib[length_reg]['GC'] == new_limit_sup)[0].tolist()))
								origin_GC_genome_approx = np.intersect1d(GC_genome,origin_genome)
							except:
								continue
#						if (origin_GC_genome.size+len(origin_GC_genome_approx))<origin_GC_pos.size+origin_GC_pos.size/5:
						origin_GC_genome_closecall=origin_pos
						#if elt == "3UTR": ## LAURA edit (before "if elt == "3UTR" and :" -- but I don't know if there was actually meant to be something after "and")
							#new_origin_genome =[i for i,v in enumerate(temp_genome_lib[length_reg]['type']) if 'intergenic' == v]
							#GC_genome = np.where(temp_genome_lib[length_reg]['GC'] == GC_searched)[0].tolist()
							##origin_GC_genome_closecall=np.intersect1d(GC_genome,new_origin_genome)
							#origin_GC_genome_closecall=new_origin_genome
						#if elt == "5UTR":
							#new_origin_genome =[i for i,v in enumerate(temp_genome_lib[length_reg]['type']) if 'promoter' == v]
							#GC_genome = np.where(temp_genome_lib[length_reg]['GC'] == GC_searched)[0].tolist()
							#origin_GC_genome_closecall=np.intersect1d(GC_genome,new_origin_genome)
							#origin_GC_genome_closecall=new_origin_genome
						if len(origin_GC_genome_approx) != 0:
							print("WARNING: not enough regions with close GC, extending to GC nearby (-"+str((i-1)*10)+"% - +"+str((i-1)*10)+"%) - "+str(elt))
							print("  this script will still search for deltaGC="+str(deltaGC)+" between a positive region and a potential negative one.")
							
					region_pos_considered = list_pos[origin_GC_pos]
					list_neg.extend(extract_neg_region(length_reg,origin_GC_genome, region_pos_considered, origin_GC_genome_approx, genome_bin_size, imin, imax,deltaGC, origin_GC_genome_closecall))
	return list_neg
################################################
'''
OUTPUT function
'''
def print_neg_regions(positions_neg_region, p):
	global output_dir
	global name_neg_file
	with open(output_dir+name_neg_file+"_"+str(p)+"_neg.bed", "w") as OUT:
#		map(lambda x: OUT.write(x+"\n"),positions_neg_region)
		for feature in positions_neg_region:
			OUT.write(feature+"\n")
	OUT.close()
################################################			 MAIN				################################################
random.seed(a = seed_rand) # seeding randomness
bin_size_list=[]
print("\n"+"="*50+"\n")
try:
	print("- Loading genome file")
	Fasta_tair = SeqIO.to_dict(SeqIO.parse(name_fasta_file, "fasta"))
	print("- Importing information file")
	info_regions = pybedtools.BedTool(name_bed_file)
	(positions_pos_regions_filename, Nmin, Nmax) = treat_positive_region(name_pos_file, name_fasta_file, name_bed_file,simpletype)
	list_pos = np.genfromtxt(positions_pos_regions_filename,dtype = [('chrom','U25'),('start',int),('stop',int),('type','U300'),('GC',float),('length',int)])
	#	list_=np.unique([v for i,v in enumerate(list_pos['type'])]).tolist()
	#	for elt in list_:
	#		print(elt,len([i for i,v in enumerate(list_pos['type']) if elt == v]))
	print("- Cutting genome into pieces - please wait")
	genome_regions_library = {}
	(genome_regions_library, bin_size_list) = cut_this_genome_in_bin(Nmin, Nmax, name_fasta_file,simpletype)
	print("- Filtering")
	(genome_regions_library,list_pos) = filter_region(genome_regions_library, list_pos, limit_types, positive_extended_temp_filename, general_temp_filename,simpletype)
	list_bin_length = genome_regions_library.keys()
	temp_genome_lib = {}
	p = 0
	while p<amount_of_negative_set: # computing indicated amount of negative sets
		p += 1
		print("- Finding negative set n°"+str(p))
		for key in list_bin_length:
			temp_genome_lib[key] = np.copy(genome_regions_library[key])
	#		temp_genome_lib[key] = np.genfromtxt(str(genome_regions_library[key]),dtype = [('chrom','S20'),('start',int),('stop',int),('type','S300'),('GC',float),('length',int)])
		temp_list_pos = np.copy(list_pos)
		positions_neg_region = fit_negative_set(temp_genome_lib, temp_list_pos, genome_bin_size,deltaGC)
		print("- Printing results")
		print("\t"+str(len(positions_neg_region))+" negative regions detected")
		print_neg_regions(positions_neg_region, p)	
	end = time.time()
	print("- Execution time: "+str(end - start)+ " secs")
	print("\n"+"="*50+"\n")
except:
	raise
finally: # removing tempfiles even if exit made in error state to prevent hard drive saturation
	try:
		os.remove(positive_temp_filename)
		os.remove(positive_temp_filename+"_2")
		os.remove(positive_extended_temp_filename)
	except:
		raise
	for elt in bin_size_list:
		os.remove(general_temp_filename+"_"+str(elt))



################################################			 END MAIN				################################################
