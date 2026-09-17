#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
compute_coverage.py
===================

DESCRIPTION:
This script calculates coverage values for genomic regions from an intersection file.
It processes a BED-like file containing intersection information and computes the
coverage area for each region by summing up the weighted contributions of overlapping segments.

The script is specifically designed to work with merged intersection data, where each line
represents a segment of a genomic region with associated coverage information.

RELATED TO pairwize_comparison

USAGE:
The script is designed to be called from the command line with the following arguments:

```bash
python compute_coverage.py --input_inter <input_file> [--merged]
```

ARGUMENTS: 
Argument,Description
- --input_inter, -i, Path to the input intersection file. Required.
- --merged, -m, Flag indicating that the input file contains merged intersection data.

OUTPUTS 
- <input_file>.cov: Contains the genomic regions with their calculated coverage values.

NOTES
- The input file should be a tab-separated file with at least 8 columns.
- Column 7 is used as the weight for coverage calculation.
- The script assumes that the input file is sorted by genomic coordinates.
- The --merged flag is currently mandatory as the script only supports merged data processing.
"""

###### LOAD ARGUMENTS ######
import argparse

parser = argparse.ArgumentParser()											   
parser.add_argument("--input_inter", "-i", type=str)
parser.add_argument("--merged", "-m",action='store_true', default= False)

args = parser.parse_args()
input_inter=args.input_inter
merged=args.merged

###### MAIN ######
with open(input_inter,"r") as f1:
	tab=[]
	for k,line in enumerate(f1): #count lines
		pass

with open(input_inter,"r") as fi:
	with open(input_inter+".cov","w") as f1:
		if not merged:
			raise Exception("option not available anymore")
		else:
			begin=[0,0,0,0]  
			area=[]
			tab=[]
			for i,line in enumerate(fi):
				if i==0:
					line=line.strip("\n").split("\t")
					tab.append(line)
					tab.append([])
					continue
				else:
					line=line.strip("\n").split("\t")
					tab[1]=list(tab[0])
					tab[0]=list(line)
					elt=list(tab[1])
				if elt[0] != begin[0] or elt[1] != begin[1] or elt[2] != begin[2] or elt[3] != begin[3]:
					begin=tab[1]
					square=(int(elt[6]) - int(elt[1]))*float(elt[7])
					if int(begin[2]) <= int(begin[6]):
						square=(int(elt[2]) - int(elt[1]))*float(elt[7])
						# area.append([elt[0],elt[1],elt[2],elt[3],square])
						f1.write(str(elt[0])+"\t"+str(elt[1])+"\t"+str(elt[2])+"\t"+str(elt[3])+"\t"+str(square)+"\n")
				else:
					if  tab[0][0] != begin[0] or tab[0][1] != begin[1] or tab[0][2] != begin[2] or  tab[0][3] != begin[3] or i==k :
						square=square+(int(elt[2]) - int(elt[5]))*float(elt[7])
						# area.append([elt[0],elt[1],elt[2],elt[3],square])
						f1.write(str(elt[0])+"\t"+str(elt[1])+"\t"+str(elt[2])+"\t"+str(elt[3])+"\t"+str(square)+"\n")
					else:
						square=square+(int(elt[6]) - int(elt[5]))*float(elt[7])

