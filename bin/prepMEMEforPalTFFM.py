#!/home/304.3-STRUCTPDEV/.local/virtualenv/MyEnv_python3/bin/python
# -*- coding: utf-8 -*-

"""
MEME Palindromic Sequence Preparer
===================================

DESCRIPTION:
This script prepares FASTA sequences for MEME analysis by processing palindromic motifs.
It reads a MEME output file, extracts sequences used to build a specific palindromic matrix,
splits each sequence in half, reverse complements the 3' half, and outputs both halves as separate sequences.
This preparation is useful for analyzing palindromic motifs where both halves of the sequence are important.

RELATED TO compute_motif

USAGE:
```bash
python script.py --meme <meme_file> --fasta <output_fasta> --pfmtoUse <pfm_index> [--seqlen <min_seq_length>]

ARGUMENTS:
- --meme, -m, Path to the MEME output file. Required.
- --fasta, -f, Path to the output FASTA file. Required.
- --pfmtoUse, -p, Index of the PFM/motif to use (0-based index). Required.
- --seqlen, -s, Minimum sequence length to keep (default: 0). Optional.

NOTES:
- The script assumes the MEME file follows the standard MEME output format.
- The script stops processing after it has read all sequences for the specified motif.
- The minimum sequence length filter (seqlen) helps exclude short sequences that might not be informative.
- The script handles sequences containing '.' characters by removing them before processing.
"""

# CONFIG
import sys
import os
import argparse
from Bio.Seq import Seq
import re

print("[WARNING] - Scripts works for full sequences")

# LOAD ARGUMENTS
parser = argparse.ArgumentParser() 
parser.add_argument("--meme","-m")
parser.add_argument("--fasta","-f")
parser.add_argument("--pfmtoUse","-p", type=int)
parser.add_argument("--seqlen","-s", type=int, default=0) 
# this depends on the meme alignment and the size of the symetric motif..
args = parser.parse_args()

# SUPLEMNETAL CONFIG 
fin=open("%s" % args.meme, "r")
fou=open("%s" % args.fasta, "w")
MINSEQLEN=int(args.seqlen)
gate="closed"
cpt=0
pfmtouse=args.pfmtoUse+1
nseq=3000

# MAIN
def micmac(LINE):
  tmp=LINE.strip("\n").split(" ")
  #print len(tmp)
  #print tmp
  while("" in tmp) : 
    tmp.remove("") 
  seq=''.join(tmp[-3:])
  s1=seq[:len(seq)]
  if "." in s1:
    print("[WARNING] - Non ATCG char in the seq:")
    # print (tmp)
  s1=s1.replace(".", "")
  if bool(re.match('^[ACTG]+$', s1))==False:
    print ("[WARNING] - Non ATCG and non . char in the seq:")
    # print (tmp)
    # print (s1)
  #s2=seq[len(seq):]
  s2rc=Seq(s1).reverse_complement()
  return [s1, s2rc]
once=False
for l in fin:
  if l[0:5]=="MOTIF":
    #print(str(l))
    #MOTIF ACACGTGT MEME-1	width =   8  sites = 386  llr = 3574  E-value = 4.8e-591
    #MOTIF CGGCGCCG MEME-3	width =   8  sites =   8  llr = 94  E-value = 2.7e+004
    #nseq=int(l.split(" ")[8])
    
    nseq=int(list(filter(None, l.split(' ')))[7])
    print('-> Line start with motif')
    
    print("%d sequences were used to built the symetric motif\n" % nseq)
    if gate=="open":
      break

  if "MEME-"+str(pfmtouse)+" sites sorted by position p-value" in l:
    print("[INFO) - Found MEME matrix to use")
    gate="open"
    seen=0
    kept=0


  if gate=="open" and cpt>3 and l=="--------------------------------------------------------------------------------":
    break
  if gate=="open" and cpt>3 and cpt<nseq+4:
    cpt+=1
    listSeq=micmac(l)
    if len(listSeq[0])>MINSEQLEN:
      fou.write(">seq%d\n%s\n" % (seen, listSeq[0]))
      kept+=1
    if len(listSeq[1])+1>MINSEQLEN:
      fou.write(">seq%drc\n%s\n" % (seen, listSeq[1]))
      kept+=1
    seen+=1
  elif gate=="open" and cpt < 4:
    cpt+=1

fou.close()
fin.close()
print ("%d seqs were seen again\n" % seen)
print ("%d seqs were kept\n" % kept)

