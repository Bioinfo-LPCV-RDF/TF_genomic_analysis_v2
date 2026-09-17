#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SeqConv Motif Discovery
=======================

DESCRIPTION:
This script uses a pre-trained SeqConv model to discover motifs in DNA sequences.
It processes input sequences, generates predictions using the model, extracts potential
motif sequences, and creates a Position Frequency Matrix (PFM) in MEME format.

SOURCE:
- https://github.com/shenwei19/SeqConv/blob/master/src/SeqConv.py

USAGE
```bash
python seqconv_motif.py <input_file> <peak_size> <name> <output_dir> <model_file>
```

ARGUMENTS:
- <input_file>   : Path to the input data file containing sequences. Required.
- <peak_size>    : Size of peaks in base pairs. Required.
- <name>         : Prefix for output files. Required.
- <output_dir>   : Directory where output files will be saved. Required.
- <model_file>   : Path to the pre-trained Keras model file. Required.

NOTES:
- The script uses one-hot encoding to represent DNA sequences.
- It extracts 24-base sequences corresponding to the highest model activations.
- The script handles both forward and reverse complement strands for motif discovery.
- The model should be a Keras model trained using the SeqConv architecture.
"""

# IMPORTS
import numpy as np
import sys
import tensorflow as tf
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import time
import os
import sys
import keras
from keras.models import load_model 
from keras import utils
from keras import layers
from keras.utils import to_categorical
from collections import Counter
from keras.models import Model
from keras.layers import Input, Dense, Dropout, Flatten
from keras.layers import Convolution2D, MaxPooling2D, AveragePooling2D, ReLU
from keras.layers import BatchNormalization
from keras.utils import plot_model
from keras.callbacks import EarlyStopping
from keras.layers import concatenate
from keras.initializers import random_normal
from keras import regularizers
from sklearn.metrics import roc_curve, auc, average_precision_score
from tqdm import tqdm

np.random.seed(12580)

# LAOD ARGUMENTS 
in_file = sys.argv[1]
peaksize = sys.argv[2]
name = sys.argv[3]
outdir = sys.argv[4]
model = sys.argv[5]
# print(in_file)
# print(outdir)
# print(name)
peaksize = int(float(peaksize))

np.random.seed(12580)

# FUNCTIONS 
# Convert sequences into matrix (one-hot encoding)
def Seq2mat(seq):
    mat1 = np.zeros([4, len(seq)])
    mat2 = np.zeros([4, len(seq)])
    for letter, index in zip(seq, range(len(seq))):
        if letter == 'A' or letter == 'a':
            mat1[0, index] = 1
            mat2[3, len(seq) - index - 1] = 1
        if letter == 'C' or letter == 'c':
            mat1[1, index] = 1
            mat2[2, len(seq) - index - 1] = 1
        if letter == 'G' or letter == 'g':
            mat1[2, index] = 1
            mat2[1, len(seq) - index - 1] = 1
        if letter == 'T' or letter == 't':
            mat1[3, index] = 1
            mat2[0, len(seq) - index - 1] = 1
        if letter == 'N':
            pass
    return mat1, mat2

# Generate a matrix of matrix
def mat2Mat(filename):
	Mat = []
	Tag = []
	Seq = []
	with open(filename) as f:
		for line in f.readlines():
			tag = line.split()[0]
			seq = line.split()[1]
			Seq.append(seq)
			mat1,mat2 = Seq2mat(seq)
			Mat.append([mat1,mat2])
			Tag=np.append(Tag,float(tag))
	return np.array(Mat),Tag,Seq

# MAIN
Mat,Tag,Seq = mat2Mat(in_file)

x_1 = Mat[:,0,:,:]
x_1 = x_1.reshape(x_1.shape[0],4,peaksize,1)
x_2 = Mat[:,1,:,:]
x_2 = x_2.reshape(x_2.shape[0],4,peaksize,1)

# MOTIF GENERATION -------------
print("- Motif generation (.meme) ...")
model = load_model(model)
model_2 = Model(inputs=model.input,outputs=[model.get_layer('re_lu').get_output_at(0),model.get_layer('re_lu_1').get_output_at(0)]) 
# Reverse complement function
def revComp(seq):
    rev_seq = ''
    for i in seq:
        if i in 'Aa':
            rev_seq = 'T' + rev_seq
        elif i in 'Cc':
            rev_seq = 'G' + rev_seq
        elif i in 'Gg':
            rev_seq = 'C' + rev_seq
        elif i in 'Tt':
            rev_seq = 'A' + rev_seq
        elif i == 'N':
            rev_seq = 'N' + rev_seq
    return rev_seq

# Sequence prediction and extraction
seq = []
for i in tqdm(range(len(x_1))):
    dic = {}
    l0 = model_2.predict([x_1[i].reshape(1, 4, peaksize, 1), x_2[i].reshape(1, 4, peaksize, 1)])[0].reshape(((peaksize-24)+1), 16)
    l1 = model_2.predict([x_1[i].reshape(1, 4, peaksize, 1), x_2[i].reshape(1, 4, peaksize, 1)])[1].reshape(((peaksize-24)+1), 16)
    if max(max(l0.flatten()), max(l1.flatten())) > 0:
        for j in range(16):
            l0_index = np.where(l0[:, j] == max(l0[:, j]))[0][0]
            l0_max = max(l0[:, j])
            l1_index = np.where(l1[:, j] == max(l1[:, j]))[0][0]
            l1_max = max(l1[:, j])

            dic['+_' + str(l0_index)] = max(dic.get('+_' + str(l0_index), 0), l0_max)
            dic['-_' + str(l1_index)] = max(dic.get('-_' + str(l1_index), 0), l1_max)

        strand, seq_index = max(dic, key=dic.get).split('_')

        if strand == '+':
            seq_24 = Seq[i][int(seq_index):int(seq_index) + 24]
        else:
            seq_24 = revComp(Seq[i])[int(seq_index):int(seq_index) + 24]
        seq.append(seq_24)

# Writing sequences to a file
with open(f'{outdir}/{name}/motif/seq.txt', 'w') as tar_file:
    for s in seq:
        tar_file.write(s + '\n')

# Generating Position Frequency Matrix (PFM)
file=open(f'{outdir}/{name}/motif/seq.txt').readlines()
tar_file=open(f'{outdir}/{name}/motif/{name}_seq.meme','w')

length=24
for i in range(length):
        list=[]
        for j in file:
                list.append(j[i])

        c = dict(Counter(list))
        sum_all = float(sum(c.values()))
        if 'A' in c.keys():
                freq_A = c['A']/sum_all
        else:
                freq_A = 0
        if 'C' in c.keys():
                freq_C = c['C']/sum_all
        else:
                freq_C = 0
        if 'G' in c.keys():
                freq_G = c['G']/sum_all
        else:
                freq_G = 0
        if 'T' in c.keys():
                freq_T = c['T']/sum_all
        else:
                freq_T = 0
        tar_file.write(str(freq_A)+'\t'+str(freq_C)+'\t'+str(freq_G)+'\t'+str(freq_T)+'\n')

tar_file.close()

