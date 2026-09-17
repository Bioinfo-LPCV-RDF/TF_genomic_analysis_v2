#!usr/bin/python

"""
FULL_METHYLATION
================

DESCRIPTION
This script is to generate a file that give the correspondance between:
  - Every single base position in a TF binding sites detected, using a PFM matrix in a TF bound region
  - The probabily of methylation of this single base
  - The binding affinity (measured as RPKM) of the peak region in DAP-seq and ampDAP-seq experiments

RELATED TO cooking_methylation

WORKFLOW
The script performs the following steps:
1. Loads chromosome sequences from a genome FASTA file.
2. Loads DAP and AMP peak RPKM data.
3. Loads methylation probability data.
4. Identifies the best scoring binding sites based on a PFM score cutoff.
5. Generates an output table with detailed information about peaks, binding sites, and methylation probabilities.

ARGUMENTS
- `-genome_File` or `-g` (str): Path to the genome FASTA file.
- `-peaksCov_File` or `-p` (str): Path to the file containing peak coverage data.
- `-methMap_File` or `-m` (str): Path to the file containing methylation probability data.
- `--scoresPFM_File` or `-s` (str): Path to the file containing PFM search results.
- `--out_dir` or `-o` (str): Directory where the output table will be saved.
- `--cutOff` or `-c` (float): PFM score threshold; only binding sites with scores above this value are considered.
- `--lengthMotif` or `-l` (int): Length of the motif in the PFM matrix.

OUTPUT
The script generates a tab-delimited output file named `table_peak_bs_meth_Zhu.txt` in the specified output directory. 
The table includes the following columns:
- Chromosome, peak start, and peak end positions.
- Coverage values for DAP and AMP experiments.
- Number of methylated sites and CG sites in the peak.
- Binding site information, including PFM score, start and end positions, strand, motif sequence, and methylation probabilities for each base in the motif.

"""

############### IMPORTS
from Bio import SeqIO
import sys
import os
import argparse


############### ARGS LOAD 
parser = argparse.ArgumentParser() 
parser.add_argument("-genome_File", "-g")   
parser.add_argument("-peaksCov_File", "-p")
parser.add_argument("-methMap_File", "-m")
parser.add_argument("--scoresPFM_File","-s")
parser.add_argument("--out_dir", "-o")
parser.add_argument("--cutOff", "-c", type=float) 
parser.add_argument("--lengthMotif", "-l", type=int) 

args = parser.parse_args()

############### CONFIG
NA="NA"
print("[INFO] - Missing value are set to %s\n" % NA)

print("\n\nPFM score cutoff is %d\n" % args.cutOff)

############### FUNCTIONS
def getNumberOfMethSitesInSegment(CHR, START, END):
  """
  Counts the number of methylation sites in a specified genomic segment where the methylation level exceeds 0.5.

  Args:
    CHR (str): Chromosome identifier.
    START (int): Start position of the genomic segment (0-based).
    END (int): End position of the genomic segment (0-based, inclusive).

  Returns:
    int: The number of methylation sites in the specified segment with methylation levels greater than 0.5.

  Notes:
    - The function assumes the existence of two global variables:
      - `set_meth`: A set containing coordinates of methylation sites in the format "CHR_POSITION".
      - `dico_meth`: A dictionary mapping coordinates to their methylation levels.
    - Coordinates are formatted as "CHR_POSITION", where CHR is the chromosome and POSITION is 1-based.
  """
  
  x=0
  for i in range(START+1, END+1):
    tmpCoords="%s_%d" % (CHR, i)
    if tmpCoords in set_meth:
      if dico_meth[tmpCoords]>0.5:
        x+=1
  return x

def getProbMethPerSite(CHR, START, END):
  """
  Calculate the probability of methylation for each site within a specified genomic range.

  Args:
    CHR (str): The chromosome identifier (e.g., "chr1", "chr2").
    START (int): The starting position (1-based index) of the genomic range.
    END (int): The ending position (1-based index) of the genomic range.

  Returns:
    list: A list of methylation probabilities (PM) for each site in the range.
        If a site is not found in the methylation data, the value NA is used.
  """
  LIST=[]
  for i in range(START,END+1):
    tmpCoords="%s_%d" % (CHR, i+1)
    if tmpCoords in set_meth:
      PM=dico_meth[tmpCoords]
    #elif dico_chr[CHR][i] in ["G", "C"]:
    #  PM=NA
    else:
      PM=NA
    LIST.append(PM)
  return LIST 


 #TODO to rm ? 
#print "Step 0: get GCs covered sites in Ecker's map per bound regions"
#Ecker="/home/304.6-RDF/LFY/LFY-pionner/data/Methylation/Ecker_lab_Nature_2013/GS#M1085222_mC_calls_Col_0.tsv"
#os.system("cat %s | grep -v chrom | awk '$6>4' | sed 's/^/chr/g' | awk '{print #$1, $2, $2}' | tr ' ' '\t' > tmpEcker" % Ecker)
#os.system("bedtools coverage -a %s -b tmpEcker | cut -f -6 > tmpCovEcker.bed" #% args.peaksCov_File)

############### PREPARE DATA AND PROCESS INPUT FILE 
# -------------
print("[STEP 1] - Loading chromosome sequences\n")
# init
dico_chr={}
f=open("%s" % args.genome_File, "r")
Cs=Gs=As=Ts=0

for r in SeqIO.parse(f, "fasta"):
  print("\t%s: %d\n" % (r.name, len(r.seq)))
  dico_chr[r.name]=r.seq
  Cs+=r.seq.count("C")
  Gs+=r.seq.count("G")
  As+=r.seq.count("A")
  Ts+=r.seq.count("T")
f.close()
#fout=open("/home/312.6-Flo_Re/312.6.1-Commun/Romain/genome_nt_count.txt", "w")
#fout.write("genome is: %s\n" % args.genome_File)
#fout.write("C: %d\n" % Cs)
#fout.write("G: %d\n" % Gs)
#fout.write("A: %d\n" % As)
#fout.write("T: %d\n" % Ts)
#fout.close()

# -------------
print("[STEP 2] - Loading DAP and AMP peak RPKM \n")

#init
dico_cov={}
f=open("%s" % args.peaksCov_File, "r")
#f=open("tmpCovEcker.bed", "r")

for l in f:
  tmp=l.strip().split("\t")
  chrom=tmp[0]
  start=int(tmp[1])
  end=int(tmp[2])

  covDAP=float(tmp[3])
  covAMP=float(tmp[4])
  #covCGs=int(tmp[5]) # this is the number of CGs that are seen in the meth map (>4 reads)
  peakCoord="%s_%d_%d" % (chrom, start, end)
  #dico_cov[peakCoord]={"covDAP":covDAP, "covAMP":covAMP, "covCGs":covCGs}
  dico_cov[peakCoord]={"covDAP":covDAP, "covAMP":covAMP}
f.close()

print("\t[REPORT] - Number of bound regions (peaks): %d\n" % len(dico_cov))

# -------------
print("[STEP 3] - Loading methylation probability\n")
dico_meth={}
f=open("%s" % args.methMap_File, "r")
for l in f:
  tmp=l.strip().split("\t")
  chrom=tmp[0]
  pos=int(tmp[1])
  probMeth=abs(float(tmp[3]))
  coords="%s_%d" % (chrom, pos)
  dico_meth[coords]=probMeth
f.close()

# -------------
print("[STEP 4] - Loading the best scoring binding sites")
print("\t[INFO] - If it has score > %d\n" % args.cutOff)
fin=open("%s" % args.scoresPFM_File, "r")

# init
dico={}
for l in fin:
  tmp=l.strip().split("\t")
  chrom=tmp[0].split(":")[0]
  peakStart=int(tmp[0].split(":")[1].split("-")[0])
  peakEnd=int(tmp[0].split(":")[1].split("-")[1])
  peak="%s_%d_%d" % (chrom, peakStart, peakEnd)
  pos=int(tmp[1])
  strand=tmp[3]
  score=float(tmp[7])
  if score > args.cutOff:
    if peak in dico.keys():
      if score > dico[peak][2]:
        dico[peak]=[pos, strand, score]
    else:
      dico[peak]=[pos, strand, score]
fin.close()
print("\t[REPORT] - There is %d peaks with at least one BS with pfm score > %d\n" % (len(dico), args.cutOff))


############### MAIN
# -------------
print("[STEP 5] - Print out the table\n")
cpt=0
setPeakWithBS=set(dico.keys())
set_meth=set(dico_meth.keys())

fout=open("%s/table_peak_bs_meth_Zhu.txt" % args.out_dir, "w")
#fout.write("chrom\tpeakStart\tpeakEnd\tcovDAP\tcovAMP\tmethSitesInPeak\tEckerCovCGsInPeak\tCGsInPeak\t")
fout.write("chrom\tpeakStart\tpeakEnd\tcovDAP\tcovAMP\tmethSitesInPeak\tCGsInPeak\t")
fout.write("BSpfmScore\tBSstart\tBSend\tstrand\tmotifSeq\tmethSitesInMotif\t")
for x in range(1,args.lengthMotif):
  fout.write("s%d\tpm%d\t" % (x, x))
fout.write("s%d\tpm%d\n" % (args.lengthMotif, args.lengthMotif))

for peak in dico_cov.keys():
  chrom=peak.split("_")[0]
  peakStart=int(peak.split("_")[1])
  peakEnd=int(peak.split("_")[2])
  methSitesInPeak=getNumberOfMethSitesInSegment(chrom, peakStart, peakEnd)
  peakSeq=dico_chr[chrom][peakStart:peakEnd]
  covDAP=dico_cov[peak]["covDAP"]
  covAMP=dico_cov[peak]["covAMP"]
  #EckerCovCGsInPeak=dico_cov[peak]["covCGs"]
  CGsInPeak=peakSeq.count("C")+peakSeq.count("G")

  #fout.write("%s\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t" % (chrom, peakStart, peakEnd, covDAP, covAMP, methSitesInPeak, EckerCovCGsInPeak, CGsInPeak))
  fout.write("%s\t%d\t%d\t%f\t%f\t%d\t%d\t" % (chrom, peakStart, peakEnd, covDAP, covAMP, methSitesInPeak, CGsInPeak))
 
  # DEAL WITH BS 
  if peak in setPeakWithBS:
    #print "BS found"
    pos=dico[peak][0]
    strand=dico[peak][1]
    BSpfmScore=dico[peak][2]
    BSstart=peakStart+pos-1
    BSend=BSstart+args.lengthMotif-1
    motifSeq=dico_chr[chrom][BSstart:BSend+1]

    motifSeqList=list(motifSeq)
    #print "\t\tget number of meth site in motif\n"
    methSitesInMotif=getNumberOfMethSitesInSegment(chrom, BSstart, BSend)
    #print "\t\tget prob of meth per base in motif\n"
    PMlist=getProbMethPerSite(chrom,BSstart, BSend)
    #print peakSeq
    #print "%s\t%d\t%d\t%f\t%f\t%d\t" % (chrom, peakStart, peakEnd, covDAP, covAMP, methSitesInPeak)
    #print "%f\t%d\t%d\t%d\t%s\t%d\t" % (BSpfmScore, BSstart, BSend, strand, motifSeq, methSitesInMotif)
    #print motifSeqList
    #print PMlist
    #print "\n\n"
    for i in range(0, args.lengthMotif):
      #print motifSeq[i]
      if motifSeq[i] in ["G", "C"] and PMlist[i]==NA:
        #exit("error") #these are missing values, not error
        pass
      elif motifSeq[i] not in ["G", "C"] and PMlist[i]!=NA:
        exit("erreur:i=%d, nt=%s, PM=%f\n" % (i, motifSeq[i], PMlist[i]))

    fout.write("%f\t%d\t%d\t%s\t%s\t%d\t" % (BSpfmScore, BSstart, BSend, strand, motifSeq, methSitesInMotif))
    for x in range(0, args.lengthMotif-1):
      fout.write("%s\t%s\t" % (motifSeqList[x], PMlist[x]))
    fout.write("%s\t%s\n" % (motifSeqList[args.lengthMotif-1], PMlist[args.lengthMotif-1]))

  else:
    pos=strand=BSpfmScore=BSstart=BSend=methSitesInMotif=motifSeq=NA
    PMlist=motifSeqList=[NA]*args.lengthMotif 
    fout.write("%s\t%s\t%s\t%s\t%s\t%s\t" % (BSpfmScore, BSstart, BSend, strand, motifSeq, methSitesInMotif))
    for x in range(0, args.lengthMotif-1):
      #print (motifSeqList[x], PMlist[x])
      fout.write("%s\t%s\t" % (motifSeqList[x], PMlist[x]))
    fout.write("%s\t%s\n" % (motifSeqList[args.lengthMotif-1], PMlist[args.lengthMotif-1]))

  cpt+=1
  if cpt%100==0:
    print("%d out of %d peaks" % (cpt, len(dico_cov)))
  #if cpt>2:
    #break
fout.close()


