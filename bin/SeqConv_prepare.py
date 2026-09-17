#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SeqConv Data Preparation Script
===============================

DESCRIPTION:
This script prepares training data for the SeqConv model by processing peak files
and generating positive and negative labeled sequences. It extracts sequences from
a genome file based on peak coordinates and creates a balanced dataset with both
positive (peak) and negative (random genomic regions) examples.

SOURCES:
- https://github.com/shenwei19/SeqConv/blob/master/src/prepare.py
- https://github.com/shenwei19/SeqConv/blob/master/src/generator.py

USAGE:
```bash
python prepare_seqconv_data.py <peaks_file> <peak_size> <genome_file> <name> <train_data_file> <size_file> <output_dir>
```

ARGUMENTS:
- <peaks_file>        : Path to the peak file (.narrowPeak or .bed format). Required.
- <peak_size>         : Size of peaks in base pairs (200 by default). Required.
- <genome_file>       : Path to the genome FASTA file. Required.
- <name>              : Prefix for output files. Required.
- <train_data_file>   : Path to the output training data file. Required.
- <size_file>         : Path to the chromosome sizes file. Required.
- <output_dir>        : Directory where output files will be saved. Required.

OUTPUTS:
- <name>/fa/<name>.fa: FASTA file containing positive sequences.
- <name>/train/<name>_train.txt: Training data file with labeled sequences.

NOTES: 
- It automatically creates necessary subdirectories (fa/, train/) in the output directory.
- Sequences are 0 or 1 tagged (positive or negative label)
- The number of negative samples generated depends on the number of positive peaks:
  - 20,020 negative samples if there are more than 3,000 positive peaks.
  - 3,000 negative samples otherwise.
- The script uses system calls to external tools (seqtk, bedtools) which must be available.
"""

# IMPORT  
import os
import sys
import numpy as np
import random

# LOAD ARGUMENTS
peaks = sys.argv[1] # *.bed - peaks file
peak_size = sys.argv[2] # 200 by default
genome = sys.argv[3] # *.fa - genome file
name = sys.argv[4] 
train_data = sys.argv[5] # *train.txt - an intermediate file
sizefile = sys.argv[6] # .chrmsize 
outdir=sys.argv[7]
neg=sys.argv[8] #negative set for training 

peak_size = int(float(peak_size))

# MAIN 
def prepare_data(infile1, infile2, name, size, chrm_size, outfile, outdir, negfile):
	
	# -> Positive tag/label data written in train.txt
	print('- Labelling of positive data')
	os.system('seqtk subseq %s %s > %s/%s/fa/%s.fa' %
			  (infile2, infile1, outdir, name, name))  # extract sequences
	
	file1 = open(f"{outdir}/{name}/fa/{name}.fa").readlines()
	with open(outfile, 'w') as tar_file:
		for i in range(int(len(file1)/2)):
			seq = file1[2*i + 1].strip()
			tar_file.write(f"1\t{seq}\n")

	# -> Negative tag/label data imported from file
	print('- Importing provided negative set')

	# Optional but safe: remove regions overlapping positive peaks
	tmp_neg = f"{outdir}/{name}/neg_cleaned.bed"
	os.system(f"bedtools intersect -a {negfile} -b {infile1} -v > {tmp_neg}")

	# Extract sequences and label them as 0
	os.system(
		f"seqtk subseq {infile2} {tmp_neg} | grep -v '>' | awk '{{print 0\"\t\"$0}}' >> {outfile}"
	)

	print('- Done! Negative set imported.')


	print('- Final formatting')
	os.system('bedtools intersect -a %s/%s/tmp.txt -b %s -v > %s/%s/tmp2.txt' % (outdir, name, infile1, outdir, name) ) #substract peaks already present in the positive set if necessary
	os.system('seqtk subseq %s %s/%s/tmp2.txt | grep -v \'>\' | awk \'{print 0"\t"$0}\' >> %s/%s/train/%s_train.txt' % (infile2, outdir, name, outdir, name, name)) #get sequencing and write peaks labelled as 0
	os.system('rm %s/%s/tmp.txt' % (outdir, name)) #clean
	
prepare_data(peaks, genome, name, peak_size, sizefile, train_data, outdir)