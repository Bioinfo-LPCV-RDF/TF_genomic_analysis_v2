#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SeqConv Model Prediction and Evaluation
========================================

DESCRIPTION:
This script converts input sequence data into matrices and uses a pre-trained SeqConv model
to predict transcription factor binding. It generates prediction scores and evaluates model
performance using ROC curves and AUC metrics.

SOURCE:
- https://github.com/shenwei19/SeqConv/blob/master/src/SeqConv.py

USAGE:
```bash
python seqconv_predict.py <input_file> <peak_size> <name> <output_dir> <model_file>
```

ARGUMENTS:
- <input_file>   : Path to the input data file containing sequences. Required.
- <peak_size>    : Size of peaks in base pairs. Required.
- <name>         : Prefix for output files. Required.
- <output_dir>   : Directory where output files will be saved. Required.
- <model_file>   : Path to the pre-trained Keras model file. Required.

OUTPUTS: 
- <name>_scores.txt: Text file with prediction scores.
- roc/<name>_scores.txt: BED file with prediction scores and genomic locations.
- roc/<name>_roc.txt: Text file with ROC curve values.
- roc/<name>_roc.png: PNG image of the ROC curve.
- roc/<name>_rec.txt: Text file with AUC and precision-recall scores.

NOTES:
- The script uses one-hot encoding to represent DNA sequences.
- It implements a dual-input architecture to process both forward and reverse complement sequences.
- The script saves prediction scores along with genomic coordinates for downstream analysis.
- The model file should be a Keras model trained using the SeqConv architecture.
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

# LOAD ARGUMENTS
in_file = sys.argv[1]
peaksize = sys.argv[2]
name = sys.argv[3]
outdir = sys.argv[4]
model = sys.argv[5]

peaksize = int(float(peaksize))

# FUNCTIONS
# SEQUENCES PREPARATION -------------
# convert sequences into 4 x matrices (one-hot-encoding)
def Seq2mat(seq):
    mat1 = np.zeros([4,len(seq)])
    mat2 = np.zeros([4,len(seq)])
    for letter,index in zip(seq,range(len(seq))):
        if letter == 'A' or letter == 'a':
            mat1[0,index] = 1
            mat2[3,len(seq)-index-1]=1
        if letter == 'C' or letter == 'c':
            mat1[1,index] = 1
            mat2[2,len(seq)-index-1]=1
        if letter == 'G' or letter == 'g':
            mat1[2,index] = 1
            mat2[1,len(seq)-index-1]=1
        if letter == 'T' or letter == 't':
            mat1[3,index] = 1
            mat2[0,len(seq)-index-1]=1
        if letter == 'N':
            pass
    return mat1,mat2

# convert train data set into matrices [tag (0 or 1) x peaksize x 4 (encoding) x n sequences]
def mat2Mat(filename):
    Mat = []
    Tag = []
    Seq = []
    Loc = []
    with open(filename) as f:
        for line in f.readlines():
            parts = line.split()
            # print(parts) #DEBUG
            tag = parts[0]
            seq = parts[1]
            chrom = parts[2]
            start = parts[3]
            end = parts[4]
            mat1, mat2 = Seq2mat(seq)
            if mat1.shape != (4, int(peaksize)) or mat2.shape != (4, int(peaksize)):
                raise ValueError(f"Incorrect shape for sequence: {seq}")
            Mat.append([mat1, mat2])
            Tag = np.append(Tag, float(tag))
            Seq.append(seq)
            Loc.append((chrom, start, end))
    return np.array(Mat), Tag, Seq, Loc

# MAIN
print("- Converting train dataset into matrices...")
Mat, Tag, Seq, Loc = mat2Mat(in_file)

# 2 matrices set are needed to feed the model (1 even lines, 1 odd lines)
x_1 = Mat[:,0,:,:]
x_1 = x_1.reshape(x_1.shape[0],4,int(peaksize),1)
x_2 = Mat[:,1,:,:]
x_2 = x_2.reshape(x_2.shape[0],4,int(peaksize),1)
print("- Data successfully converted")

# MODEL PREDICTION -------------
print("- Model predictions...")
model = load_model(model)
# model.summary()     
prob = model.predict([x_1, x_2])

# scores and associated ROC generation 
print("- ROC generation")
fpr, tpr, threshold = roc_curve(Tag, prob)
roc_auc = auc(fpr, tpr)
prc = average_precision_score(Tag, prob)
os.system('echo %s %s %s >> %s/roc/%s_rec.txt' % (name,roc_auc,prc,outdir,name))

np.savetxt(f'{outdir}/roc/{name}_scores.txt', prob, delimiter='\t')

with open(f'{outdir}/roc/{name}_scores.txt', 'w') as bed_file:
    for tag, loc, prob_val in zip(Tag, Loc, prob):
        chrom, start, end = loc
        bed_file.write(f'{tag}\t{chrom}\t{start}\t{end}\t.\t{prob_val[0]}\t.\n')

roc_file = open('%s/roc/%s_roc.txt' % (outdir, name),'w')
roc_file.write('fpr\ttpr\n')
for i,j in zip(fpr,tpr):
    roc_file.write(str(i)+'\t'+str(j)+'\n')
roc_file.close()

plt.figure()
plt.plot(fpr,tpr,label='area = %s' % roc_auc)
plt.title('ROC')
plt.savefig("%s/roc/%s_roc.png" % (outdir, name))

