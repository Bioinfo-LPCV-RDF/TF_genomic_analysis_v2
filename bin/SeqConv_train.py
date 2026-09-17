#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SeqConv Model Training and Evaluation
======================================

DESCRIPTION:
This script converts training data into matrices for the SeqConv model, trains the model,
and evaluates its performance. It implements a convolutional neural network for sequence
classification, specifically designed for genomic sequence analysis.

SOURCE:
- https://github.com/shenwei19/SeqConv/blob/master/src/SeqConv.py

USAGE:
```bash
python seqconv_train.py <input_file> <peak_size> <name> <output_dir>
```

ARGUMENTS:
- <input_file>    : Path to the training data file. Required.
- <peak_size>     : Size of peaks in base pairs. Required.
- <name>          : Prefix for output files. Required.
- <output_dir>    : Directory where output files will be saved. Required.

OUTPUTS: 
- <name>/loss/<name>_loss.png: Plot of training and validation loss.
- <name>/loss/<name>_loss.txt: Text file with loss values.
- <name>/roc/<name>_roc.png: ROC curve visualization.
- <name>/roc/<name>_roc.txt: Text file with ROC values.
- <name>/roc/<name>_scores.bed: BED file with prediction scores and locations.
- <name>/roc/<name>_rec.txt: Text file with AUC and precision-recall scores.
- <name>/model/<name>_model.keras: Trained Keras model file.

NOTES: 
- The script uses one-hot encoding to represent DNA sequences.
- It implements a dual-input architecture to process both forward and reverse complement sequences.
- Early stopping is used to prevent overfitting during training.
"""

# IMPORTS 
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import time
import os
import sys
import keras
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

np.random.seed(12580)

# LOAD ARGUMENTS
in_file = sys.argv[1]
peaksize = sys.argv[2]
name = sys.argv[3]
outdir = sys.argv[4]
peaksize = int(float(peaksize))

# FUNCTIONS
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
            parts = line.strip().split('\t')
            tag = parts[0]
            seq = parts[1]
            # Pad the sequence with "N" if it is shorter than peaksize
            if len(seq) < peaksize:
                gap = peaksize - len(seq)
                seq = seq + "N" * gap
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

# Split dataset into a training set and testing set randomly 
x_train, x_test, y_train, y_test, seq_train, seq_test, loc_train, loc_test = train_test_split(Mat, Tag, Seq, Loc, test_size=0.1, random_state=12580)

# 2 matrices set are needed to feed the model (1 even lines, 1 odd lines)
x_train_1 = x_train[:,0,:,:]
x_train_1 = x_train_1.reshape(x_train_1.shape[0],4,int(peaksize),1)
x_train_2 = x_train[:,1,:,:]
x_train_2 = x_train_2.reshape(x_train_2.shape[0],4,int(peaksize),1)

x_test_1 = x_test[:,0,:,:]
x_test_1 = x_test_1.reshape(x_test_1.shape[0],4,int(peaksize),1)
x_test_2 = x_test[:,1,:,:]
x_test_2 = x_test_2.reshape(x_test_2.shape[0],4,int(peaksize),1)

mat_1 = Input(shape=(4,int(peaksize),1))
mat_2 = Input(shape=(4,int(peaksize),1))
print("- Data successfully converted")

# MODEL TRAINING: ---------------------------------------------------------------------------------------
print("- Start to train the model...") 
shared_conv = Convolution2D(filters=16,kernel_size=(4,24),padding='valid') #no need to claim input_shape?

x1 = shared_conv(mat_1)
x2 = shared_conv(mat_2)

x1 = ReLU()(x1)
x2 = ReLU()(x2)

x1_1 = MaxPooling2D(pool_size=(1,int(peaksize)-24))(x1) #length-filter_length+1
x1_2 = AveragePooling2D(pool_size=(1,int(peaksize)-24))(x1)

x2_1 = MaxPooling2D(pool_size=(1,int(peaksize)-24))(x2)
x2_2 = AveragePooling2D(pool_size=(1,int(peaksize)-24))(x2)

merged_vector = concatenate([x1_1,x1_2,x2_1,x2_2])

x = Flatten()(merged_vector)
x = BatchNormalization()(x)
x = Dense(64,activation='relu',kernel_initializer=random_normal(mean=0,stddev=1),bias_initializer=random_normal(mean=0,stddev=1))(x)
main_output = Dense(1,activation='sigmoid',kernel_initializer=random_normal(mean=0,stddev=1),bias_initializer=random_normal(mean=0,stddev=1))(x)

model = Model(inputs=[mat_1,mat_2],outputs=main_output)
model.compile(loss='mse',optimizer='sgd')
early_stopping = EarlyStopping(monitor='val_loss',patience=100)

start = time.time()

history = model.fit([x_train_1,x_train_2],y_train,batch_size=32,epochs=1000,validation_data=[[x_test_1,x_test_2],y_test],callbacks=[early_stopping])

end = time.time()
dur = end - start
os.system('echo %s %s >> time_rec.txt' % (name,dur))

# --------------------------------------------------------------------------------------------------------

print("- Model testing...")

# loss computation and plot 
loss = history.history['loss']
val_loss = history.history['val_loss']

plt.figure()
plt.subplots_adjust(wspace=0.5,hspace=0.5)
plt.subplot(1,2,1)
plt.plot(range(1,len(loss)+1),loss,color='red')
plt.title('loss',fontsize=10)

plt.subplot(1,2,2)
plt.plot(range(1,len(val_loss)+1),val_loss,color='green')
plt.title('val_loss',fontsize=10)
plt.savefig("%s/%s/loss/%s_loss.png" % (outdir, name, name))

loss_file = open('%s/%s/loss/%s_loss.txt' % (outdir, name, name),'w')
loss_file.write('loss\tval_loss\n')
for i,j in zip(loss,val_loss):
    loss_file.write(str(i)+'\t'+str(j)+'\n')
loss_file.close()
print("- Loss values done!")

# ROC computation 
prob = model.predict([x_test_1,x_test_2])
fpr, tpr, threshold = roc_curve(y_test,prob)
roc_auc = auc(fpr,tpr)
prc = average_precision_score(y_test,prob)
os.system('echo %s %s %s >> %s/%s/roc/%s_rec.txt' % (name,roc_auc,prc,outdir,name,name))

i=0 #writting out scores values and their associated original peaks and label
with open(f'{outdir}/{name}/roc/{name}_scores.bed', 'w') as bed_file:
    for tag, loc, prob in zip(y_test, loc_test, prob):
        chrom, start, end = loc
        bed_file.write(f'{y_test[i]}\t{chrom}\t{start}\t{end}\t.\t{prob[0]}\t.\n')
        i+=1

roc_file = open('%s/%s/roc/%s_roc.txt' % (outdir, name, name),'w')
roc_file.write('fpr\ttpr\n')
for i,j in zip(fpr,tpr):
    roc_file.write(str(i)+'\t'+str(j)+'\n')
roc_file.close()

plt.figure()
plt.plot(fpr,tpr,label='area = %s' % roc_auc)
plt.title('ROC')
plt.savefig("%s/%s/roc/%s_roc.png" % (outdir, name, name))
print("- ROC done!")

# saving the model
model.save('%s/%s/model/%s_model.keras' % (outdir, name, name))
