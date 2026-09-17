

# -*- coding: utf-8 -*-
"""
merge_all_peaks_vextended.py
============================

DESCRIPTION
This script is designed to merge genomic peaks from two input files based on specific criteria. 
It processes the peaks to identify overlaps, containment, and proximity, and generates output files 
with unique peaks, merged peaks, and processed peaks in various formats.

RELATED TO pairwize_comparison

USAGE
	python merge_all_peaks_vextended.py -f1 <file1> -f2 <file2> -o <output_directory> [-d <data_directory>]
Arguments:
	-f1, --file1       Path to the first input file containing genomic peaks.
	-f2, --file2       Path to the second input file containing genomic peaks.
	-o, --output       Path to the output directory where results will be saved.
	-d, --data         (Optional) Path to the data directory. If not provided, it will be derived 
					   from the output directory and input file names.
OUTPUTS
	1. <data_directory>/<file1>_<file2>_peaks_uniques.bed:
	   Contains unique peaks that were not merged with any other peak.
	2. <data_directory>/<file1>_<file2>_peaks_merged.bed:
	   Contains merged peaks based on the defined criteria.
	3. <data_directory>/<file1>_<file2>_peaks_processed.gff:
	   Contains processed peaks in GFF format with detailed annotations.

FEATURES
	- Reads input peak files and stores them in a structured list.
	- Implements a dual-loop sliding window approach to compare peaks.
	- Merges peaks based on the following criteria:
		- Proximity of maximum values within a threshold distance.
		- Containment of one peak within another.
		- Partial overlap between peaks.
	- Outputs merged and unique peaks in BED and GFF formats.
	- Provides detailed annotations for merged peaks in the GFF output.

FUNCTIONS
	- check_merge(maxI, closemax, candidateclosemax):
		Determines whether a candidate peak is closer to the current peak than the existing closest peak.

PROCESS
	- Peaks are compared in a pairwise manner.
	- Merging decisions are made based on overlap, containment, and proximity criteria.
	- Peaks are annotated with reasons for merging or marked as unique if not merged.

NOTES
	- The script assumes that the input files are in BED format with specific columns for chromosome, 
	  start, stop, and other attributes.
	- The output directory must be writable, and sufficient disk space should be available for the results.
"""

print("")
print("====> Launching merge_all_peaks_vextended.py from bin/" )
print("")

from os.path import basename
from os.path import join
import argparse
from argparse import RawTextHelpFormatter

####### LOAD ARGUMENTS #######
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f1","--file1", type=str, help='file 1')
parser.add_argument("-f2","--file2", type=str, help='file 2')
parser.add_argument("-o","--output", type=str, help='output directory')
parser.add_argument("-d","--data", type=str, default="", help='data directory')

args = parser.parse_args()

file1 = args.file1
file2 = args.file2
output = args.output
data=args.data

######## CONFIG DATA DIR ########
if data == "":
	data = output+"/"+file1+"_"+file2+"/"
else:
	data=data+"/"
print("[INFO] - Data directory to be used: "+data)

######## LOAD DATA ########
"""
Reading peaks file and storing values in list
"""
all_peaks=[]
with open(data+file1+"_"+file2+"_peaks.bed","r") as f1: # open merged peaks file
	for line in f1:
		tmpline={}
		line=line.strip("\n")
		l=line.split("\t")
		# print(str(l)) # DEBUG
		tmpline["chrom"]=l[0]
		tmpline["startext"]=int(l[1])
		tmpline["stopext"]=int(l[2])
		tmpline["origin"]=l[3]
		tmpline["start"]=int(l[4])
		tmpline["stop"]=int(l[5])

		tmpline["max"]=int(l[4])+int((tmpline["stop"]-tmpline["start"])/2)
		tmpline["merged"]=False
		tmpline["relatedpeak"]=-1
		tmpline["reason"]="None"
		tmpline["passed"]=False
		all_peaks.append(tmpline)

######## FUNCTIONS ########

def check_merge(maxI,closemax,candidateclosemax):
	"""
	Determines whether the candidate close maximum is a better match 
	to the given maximum index compared to the current close maximum.

	Args:
		maxI (int): The index of the maximum value.
		closemax (int): The index of the current close maximum.
		candidateclosemax (int): The index of the candidate close maximum.

	Returns:
		bool: True if the candidate close maximum is closer to the 
			  maximum index than the current close maximum, False otherwise.
		"""
	oldprefdist=abs(maxI-closemax)
	actualprefdist=abs(maxI-candidateclosemax)
	if oldprefdist > actualprefdist:
		return True
	else:
		return False
# chr3:17,390,929..17,396,584

######## MAIN SCRIPT ########
# loop initiation 
i=0 #current line 
j=0 #next 

# parse all peaks
while (i < len(all_peaks) -1):
	j=i+1 #
	# print(str(all_peaks[i]))
	while (j < len(all_peaks)): # dual loop for sliding on the lines / start J >= start I
		# print(str(all_peaks[j]))
		# print(all_peaks[i]["chrom"],all_peaks[j]["chrom"],all_peaks[i]["merged"], all_peaks[i]["origin"], all_peaks[j]["origin"])
		if all_peaks[i]["chrom"] == all_peaks[j]["chrom"] and all_peaks[i]["origin"] != all_peaks[j]["origin"] and all_peaks[j]["startext"] < all_peaks[i]["stopext"]:
			thresholdDistanceMaximums=min(abs(all_peaks[j]["stop"]-all_peaks[i]["stop"]),abs(all_peaks[j]["start"]-all_peaks[i]["start"]))
			
			# print(str(thresholdDistanceMaximums))
			if abs(all_peaks[i]["max"]-all_peaks[j]["max"]) <= thresholdDistanceMaximums and thresholdDistanceMaximums < min(all_peaks[i]["stop"]-all_peaks[i]["start"],all_peaks[j]["stop"]-all_peaks[j]["start"])*0.25 : # maximum close to each other
				# print("\tmaxclose")
				if all_peaks[i]["merged"] is False and all_peaks[j]["merged"] is False:
					all_peaks[i]["merged"]=True
					all_peaks[j]["merged"]=True
					all_peaks[i]["reason"]="maxClose"
					all_peaks[j]["reason"]="maxClose"
				else:
					Imustmerge=True
					Jmustmerge=True
					if all_peaks[i]["merged"] is True:
						Imustmerge=check_merge(all_peaks[i]["max"],all_peaks[all_peaks[i]["relatedpeak"]]["max"],all_peaks[j]["max"])
					if all_peaks[j]["merged"] is True:
						Jmustmerge=check_merge(all_peaks[j]["max"],all_peaks[all_peaks[j]["relatedpeak"]]["max"],all_peaks[j]["max"])
					if Imustmerge and Jmustmerge:
						all_peaks[i]["merged"]=True
						all_peaks[j]["merged"]=True
						all_peaks[i]["reason"]="maxClose"
						all_peaks[j]["reason"]="maxClose"
						# resetting old merge peaks
						if all_peaks[i]["relatedpeak"] != -1:
							all_peaks[all_peaks[i]["relatedpeak"]]["merged"]=False
							all_peaks[all_peaks[i]["relatedpeak"]]["reason"]="None"
							all_peaks[all_peaks[i]["relatedpeak"]]["relatedpeak"]=-1
						else:
							all_peaks[all_peaks[j]["relatedpeak"]]["merged"]=False
							all_peaks[all_peaks[j]["relatedpeak"]]["reason"]="None"
							all_peaks[all_peaks[j]["relatedpeak"]]["relatedpeak"]=-1
					
			elif all_peaks[j]["stopext"] <= all_peaks[i]["stopext"]: # j contained in i
				# print("\tcontain")
				if all_peaks[i]["merged"] is False and all_peaks[j]["merged"] is False:
					all_peaks[i]["merged"]=True
					all_peaks[j]["merged"]=True
					all_peaks[i]["reason"]="contain"
					all_peaks[j]["reason"]="contained"
				else:
					Imustmerge=False
					Jmustmerge=False
					if all_peaks[i]["merged"] is True:
						Imustmerge=check_merge(all_peaks[i]["max"],all_peaks[all_peaks[i]["relatedpeak"]]["max"],all_peaks[j]["max"])
					if all_peaks[j]["merged"] is True:
						Jmustmerge=check_merge(all_peaks[j]["max"],all_peaks[all_peaks[j]["relatedpeak"]]["max"],all_peaks[j]["max"])
					if Imustmerge and Jmustmerge:
						all_peaks[i]["merged"]=True
						all_peaks[j]["merged"]=True
						all_peaks[i]["reason"]="maxClose"
						all_peaks[j]["reason"]="maxClose"
						# resetting old merge peaks
						if all_peaks[i]["relatedpeak"] != -1:
							all_peaks[all_peaks[i]["relatedpeak"]]["merged"]=False
							all_peaks[all_peaks[i]["relatedpeak"]]["reason"]="None"
							all_peaks[all_peaks[i]["relatedpeak"]]["relatedpeak"]=-1
						else:
							all_peaks[all_peaks[j]["relatedpeak"]]["merged"]=False
							all_peaks[all_peaks[j]["relatedpeak"]]["reason"]="None"
							all_peaks[all_peaks[j]["relatedpeak"]]["relatedpeak"]=-1



			elif all_peaks[j]["startext"] <  all_peaks[i]["stopext"] and all_peaks[j]["stopext"] > all_peaks[i]["stopext"]: # I and J overlaps each other, at least partially
				minilength=min(all_peaks[i]["stopext"] - all_peaks[i]["startext"], all_peaks[j]["stopext"] - all_peaks[j]["startext"])
				if (all_peaks[i]["stopext"] - all_peaks[j]["startext"]) /float(minilength) >= 0.75: # overlaps is >= 75% of the smaller peak
					# print("\tpartial")
					if all_peaks[i]["merged"] is False and all_peaks[j]["merged"] is False:
						all_peaks[i]["merged"]=True
						all_peaks[j]["merged"]=True
						all_peaks[i]["reason"]="partial"
						all_peaks[j]["reason"]="partial"
					else:
						Imustmerge=False
						Jmustmerge=False
						if all_peaks[i]["merged"] is True:
							Imustmerge=check_merge(all_peaks[i]["max"],all_peaks[all_peaks[i]["relatedpeak"]]["max"],all_peaks[j]["max"])
						if all_peaks[j]["merged"] is True:
							Jmustmerge=check_merge(all_peaks[j]["max"],all_peaks[all_peaks[j]["relatedpeak"]]["max"],all_peaks[j]["max"])
						if Imustmerge and Jmustmerge:
							all_peaks[i]["merged"]=True
							all_peaks[j]["merged"]=True
							all_peaks[i]["reason"]="maxClose"
							all_peaks[j]["reason"]="maxClose"
							# resetting old merge peaks
							if all_peaks[i]["relatedpeak"] != -1:
								all_peaks[all_peaks[i]["relatedpeak"]]["merged"]=False
								all_peaks[all_peaks[i]["relatedpeak"]]["reason"]="None"
								all_peaks[all_peaks[i]["relatedpeak"]]["relatedpeak"]=-1
							else:
								all_peaks[all_peaks[j]["relatedpeak"]]["merged"]=False
								all_peaks[all_peaks[j]["relatedpeak"]]["reason"]="None"
								all_peaks[all_peaks[j]["relatedpeak"]]["relatedpeak"]=-1
			all_peaks[i]["relatedpeak"]=j
			all_peaks[j]["relatedpeak"]=i
		elif all_peaks[i]["chrom"] != all_peaks[i]["chrom"] or all_peaks[i]["merged"] is True or all_peaks[j]["startext"] > all_peaks[i]["stopext"]:
			break
		j+=1
	i+=1
	if  i%1000 == 0 :
		 print(i)

with open(data+file1+"_"+file2+"_peaks_uniques.bed","w") as f1:
	for elt in all_peaks:
		if elt["merged"] is False:
			# f1.write(str(elt)+"\n")
			f1.write(str(elt["chrom"])+"\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t"+str(elt["origin"])+"\t"+str(elt["max"])+"\tNone\n")

with open(data+file1+"_"+file2+"_peaks_merged.bed","w") as f1:
	for elt in all_peaks: # pour chaque peak
		if elt["merged"] is True:
			elt2=all_peaks[elt["relatedpeak"]]
			newmax=int((elt["max"]+elt2["max"])/2)
			if elt["reason"]== "contain":
				f1.write(str(elt["chrom"])+"\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\tboth\t"+str(newmax)+"\t"+elt["reason"]+"\n")
				elt["passed"] = True
				elt2["passed"] = True
			elif elt["reason"]== "contained":
				f1.write(str(elt["chrom"])+"\t"+str(elt2["startext"])+"\t"+str(elt2["stopext"])+"\tboth\t"+str(newmax)+"\t"+elt2["reason"]+"\n")
				elt["passed"] = True
				elt2["passed"] = True
			elif elt["passed"] != True and elt2["passed"] != True:
				start=min(elt["startext"],elt2["startext"])
				stop=max(elt["stopext"],elt2["stopext"])
				f1.write(str(elt["chrom"])+"\t"+str(start)+"\t"+str(stop)+"\tboth\t"+str(newmax)+"\t"+elt["reason"]+"\n")
				elt["passed"] = True
				elt2["passed"] = True
# quit()
with open(data+file1+"_"+file2+"_peaks_processed.gff","w") as f1:
		IDpeaksCounter = 1
		for elt in all_peaks: # pour chaque peak
			if elt["merged"] is True:
				elt2=all_peaks[elt["relatedpeak"]]
				newmax=int((elt["max"]+elt2["max"])/2)
				if elt["reason"]== "contain":
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tboundRegion\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+";full_name=both\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tPeak\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".1;Parent=both."+str(IDpeaksCounter)+"\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tMaximum\t"+str(elt["max"]-1)+"\t"+str(elt["max"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".1.1;Parent=both."+str(IDpeaksCounter)+".1\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tPeak\t"+str(elt2["startext"])+"\t"+str(elt2["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".2;Parent=both."+str(IDpeaksCounter)+"\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tMaximum\t"+str(elt2["max"]-1)+"\t"+str(elt2["max"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".2.1;Parent=both."+str(IDpeaksCounter)+".2\n")
				elif elt["reason"]== "contained":
					f1.write(str(elt["chrom"])+"\t"+str(elt2["reason"])+"\tboundRegion\t"+str(elt2["startext"])+"\t"+str(elt2["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+";full_name=both\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt2["reason"])+"\tPeak\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".1;Parent=both."+str(IDpeaksCounter)+"\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt2["reason"])+"\tMaximum\t"+str(elt["max"]-1)+"\t"+str(elt["max"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".1.1;Parent=both."+str(IDpeaksCounter)+".1\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt2["reason"])+"\tPeak\t"+str(elt2["startext"])+"\t"+str(elt2["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".2;Parent=both."+str(IDpeaksCounter)+"\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt2["reason"])+"\tMaximum\t"+str(elt2["max"]-1)+"\t"+str(elt2["max"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".2.1;Parent=both."+str(IDpeaksCounter)+".2\n")
				elif elt["reason"] != "None" and elt2["reason"] != "None":
					start=min(elt["startext"],elt2["startext"])
					stop=max(elt["stopext"],elt2["stopext"])
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tboundRegion\t"+str(start)+"\t"+str(stop)+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+";full_name=both\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tPeak\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".1;Parent=both."+str(IDpeaksCounter)+"\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tMaximum\t"+str(elt["max"]-1)+"\t"+str(elt["max"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".1.1;Parent=both."+str(IDpeaksCounter)+".1\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tPeak\t"+str(elt2["startext"])+"\t"+str(elt2["stopext"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".2;Parent=both."+str(IDpeaksCounter)+"\n")
					f1.write(str(elt["chrom"])+"\t"+str(elt["reason"])+"\tMaximum\t"+str(elt2["max"]-1)+"\t"+str(elt2["max"])+"\t.\t.\t.\tID=both."+str(IDpeaksCounter)+".2.1;Parent=both."+str(IDpeaksCounter)+".2\n")
				elt["reason"] = "None"
				elt2["reason"] = "None"
			else:
				f1.write(str(elt["chrom"])+"\tNotMerged\tboundRegion\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t.\t.\t.\tID="+str(elt["origin"])+"."+str(IDpeaksCounter)+";full_name="+str(elt["origin"])+"\n")
				f1.write(str(elt["chrom"])+"\tNotMerged\tPeak\t"+str(elt["startext"])+"\t"+str(elt["stopext"])+"\t.\t.\t.\tID="+str(elt["origin"])+"."+str(IDpeaksCounter)+".1;Parent="+str(elt["origin"])+"."+str(IDpeaksCounter)+"\n")
				f1.write(str(elt["chrom"])+"\tNotMerged\tMaximum\t"+str(elt["max"]-1)+"\t"+str(elt["max"])+"\t.\t.\t.\tID="+str(elt["origin"])+"."+str(IDpeaksCounter)+".1.1;Parent="+str(elt["origin"])+"."+str(IDpeaksCounter)+".1\n")
			IDpeaksCounter+=1


# Chr01	.	gene	382865	384790	.	-	.	ID=Cs01G00001
# Chr01	.	mRNA	382865	384790	.	-	.	ID=Cs01G00001.1;Parent=Cs01G00001
# Chr01	.	five_prime_UTR	384657	384790	.	-	.	ID=Cs01G00001.1.utr5p1;Parent=Cs01G00001.1
# Chr01	.	exon	384339	384790	.	-	.	ID=Cs01G00001.1.exon1;Parent=Cs01G00001.1
# Chr01	.	CDS	384339	384656	.	-	0	ID=cds.Cs01G00001.1;Parent=Cs01G00001.1
# Chr01	.	exon	383804	383922	.	-	.	ID=Cs01G00001.1.exon2;Parent=Cs01G00001.1
# Chr01	.	CDS	383804	383922	.	-	0	ID=cds.Cs01G00001.1;Parent=Cs01G00001.1
# Chr01	.	exon	382865	383681	.	-	.	ID=Cs01G00001.1.exon3;Parent=Cs01G00001.1
# Chr01	.	CDS	383180	383681	.	-	1	ID=cds.Cs01G00001.1;Parent=Cs01G00001.1
# Chr01	.	three_prime_UTR	382865	383179	.	-	.	ID=Cs01G00001.1.utr3p1;Parent=Cs01G00001.1
# Chr01	.	gene	386877	389683	.	+	.	ID=Cs01G00002
# Chr01	.	mRNA	386877	389683	.	+	.	ID=Cs01G00002.1;Parent=Cs01G00002
# Chr01	.	five_prime_UTR	386877	387275	.	+	.	ID=Cs01G00002.1.utr5p1;Parent=Cs01G00002.1
# Chr01	.	exon	386877	387387	.	+	.	ID=Cs01G00002.1.exon1;Parent=Cs01G00002.1
# Chr01	.	CDS	387276	387387	.	+	0	ID=cds.Cs01G00002.1;Parent=Cs01G00002.1
# Chr01	.	exon	387552	387733	.	+	.	ID=Cs01G00002.1.exon2;Parent=Cs01G00002.1
# Chr01	.	CDS	387552	387733	.	+	2	ID=cds.Cs01G00002.1;Parent=Cs01G00002.1
# Chr01	.	exon	387891	388080	.	+	.	ID=Cs01G00002.1.exon3;Parent=Cs01G00002.1
# Chr01	.	CDS	387891	388080	.	+	0	ID=cds.Cs01G00002.1;Parent=Cs01G00002.1
# Chr01	.	exon	388157	388548	.	+	.	ID=Cs01G00002.1.exon4;Parent=Cs01G00002.1
# Chr01	.	CDS	388157	388548	.	+	2	ID=cds.Cs01G00002.1;Parent=Cs01G00002.1
# Chr01	.	exon	388635	389683	.	+	.	ID=Cs01G00002.1.exon5;Parent=Cs01G00002.1
# Chr01	.	CDS	388635	389411	.	+	0	ID=cds.Cs01G00002.1;Parent=Cs01G00002.1
# Chr01	.	three_prime_UTR	389412	389683	.	+	.	ID=Cs01G00002.1.utr3p1;Parent=Cs01G00002.1


# chr1	Araport11	gene	3631	5899	.	+	.	ID=AT1G01010;Name=AT1G01010;full_name=NAC domain containing protein 1;computational_description=NAC domain containing protein 1;locus=2200935;symbol=NAC001;Note=NAC domain containing protein 1;Dbxref=TAIR:2200935;locus_type=protein_coding
# chr1	Araport11	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;Name=AT1G01010.1;Parent=AT1G01010;full_name=NAC domain containing protein 1;computational_description=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;symbol=NAC001;Note=NAC domain containing protein 1;conf_rating=****;Dbxref=TAIR:2200934,UniProt:Q0WV96
# chr1	Araport11	CDS	3760	3913	.	+	0	ID=AT1G01010:CDS:1;Name=NAC001:CDS:1;Parent=AT1G01010.1;computational_description=NAC domain containing protein 1;curator_summary=Member of the NAC domain containing family of plant specific transcriptional regulators.;Note=NAC domain containing protein 1

quit()

# -*- coding: utf-8 -*-

from os.path import basename
from os.path import join
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f1","--file1", type=str, help='file 1')
parser.add_argument("-f2","--file2", type=str, help='file 2')
parser.add_argument("-o","--output", type=str, help='output directory')
parser.add_argument("-d","--data", type=str, default="", help='data directory')

args = parser.parse_args()

file1 = args.file1
file2 = args.file2
output = args.output
data=args.data

# data directory prep
if data == "":
	data = output+"/"+file1+"_"+file2+"/"
else:
	data=data+"/"
print("Data directory to be used: "+data)

"""
Reading peaks file and storing values in list
"""
all_peaks=[]
with open(data+file1+"_"+file2+"_peaks.bed","r") as f1:
	for line in f1:
		line=line.strip("\n")
		l=line.split("\t")
		# print(line)
		l[1]=int(l[1]) # start
		l[2]=int(l[2]) # end
		l[4]=int(l[4]) # start
		l[5]=int(l[5]) # end
		l.append(False) # = l[6]
		all_peaks.append(l) # l -> chr start end name False
		
i=0
j=0
while (i < len(all_peaks) -1):
	j=i+1
	# print(str(all_peaks[i]))
	while (j < len(all_peaks)): # dual loop for sliding on the lines / start J >= start I
		# print("\t"+str(all_peaks[j]))
		# Start J >= start I
		if all_peaks[i][0] == all_peaks[j][0] and all_peaks[i][6] is False and all_peaks[i][3]!=all_peaks[j][3]: # même CHR # TODO enlever condition True ?
			# print("\t\tstep1")
			print(str(all_peaks[j][2])+" "+str(all_peaks[i][2]))
			if  all_peaks[j][2] <= all_peaks[i][2] : # stop de J avant stop I / I contient J
				# print("\t\tJ in i")
				all_peaks[j][6]=True # a.k.a. merged
				all_peaks[i][6]=True # a.k.a. merged
				all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3],all_peaks[j][1],all_peaks[j][2],1]) 
				if (all_peaks[j][1] - all_peaks[i][1]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 and len(all_peaks[i][7])<2: 
					# If length of I before start of J higher than 50% of the peak -> create a new spe peaks out of this - not re tested for merge !
					all_peaks.insert(i,[all_peaks[i][0],all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],all_peaks[i][1],all_peaks[j][1]-1,False])
					# all_peaks.insert(i,["bis",all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],False])
					i+=1
					j+=1
				if (all_peaks[i][2] - all_peaks[j][2]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5  and len(all_peaks[i][7])<2:
					# If length of I after end of J higher than 50% of the peak -> create a new spe peaks out of this  - re tested for potential merge
					all_peaks.insert(j+1,[all_peaks[i][0],all_peaks[j][2]+1,all_peaks[i][2],all_peaks[i][3],all_peaks[j][2]+1,all_peaks[i][2],False])

			elif ((all_peaks[j][1] < all_peaks[i][2]) and (all_peaks[j][2] > all_peaks[i][2])) : # start de J contenu dans I mais pas le stop.
				# print("\t\tpartial overlap")
				if all_peaks[i][2] -all_peaks[i][1] < all_peaks[j][2] -all_peaks[j][1]: # I plus petit que J
					if (all_peaks[i][2] - all_peaks[j][1]) /float(all_peaks[i][2] - all_peaks[i][1]) > 0.66:
						all_peaks[j][6]=True # a.k.a. merged
						all_peaks[i][6]=True # a.k.a. merged
						all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3],all_peaks[j][1],all_peaks[j][2],2]) 
						if (all_peaks[j][2] - all_peaks[i][2]) / float(all_peaks[j][2] - all_peaks[j][1] ) > 0.5  and len(all_peaks[i][7])<2: # Non overlap de J / taille de J > 0.5 
							all_peaks.insert(j+1,[all_peaks[i][0],all_peaks[i][2]+1,all_peaks[j][2],all_peaks[j][3],all_peaks[i][2]+1,all_peaks[j][2],False])
				else: # I plus grand ou égal à J
					if (all_peaks[i][2] - all_peaks[j][1]) /float(all_peaks[j][2] - all_peaks[j][1]) > 0.66: # overlap / taille de J > 0.7
						all_peaks[j][6]=True # a.k.a. merged
						all_peaks[i][6]=True # a.k.a. merged
						all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3],all_peaks[j][1],all_peaks[j][2],3])
						if (all_peaks[j][1] - all_peaks[i][1]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5  and len(all_peaks[i][7])<2: # Non overlap de I / taille de I > 0.5
							all_peaks.insert(i,[all_peaks[i][0],all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],all_peaks[i][1],all_peaks[j][1]-1,False])
							i+=1
							j+=1						
			elif ((all_peaks[j][1] >= all_peaks[i][2])): #JS après IE -> pas de merge, on peut arreter le while il n'y aura plus de peak à merge après
				break
			if j+1<len(all_peaks):
				if ((all_peaks[i][2] >= all_peaks[j+1][1]) and all_peaks[i][3] != all_peaks[j+1][3]):
					all_peaks[i][6]=False
		elif all_peaks[i][3]!=all_peaks[j][3] : # CHR différents, plus de merge possible après
			break
		# if j-i == 2:
		#	 print j-i,all_peaks[i][0],all_peaks[i][1],all_peaks[i][2],all_peaks[j][1],all_peaks[j][2]
		j+=1
	i+=1
	if  i%1000 ==0 :
		 print(i)
		
with open(data+file1+"_"+file2+"_peaks_uniques.bed","w") as f1:
	for elt in all_peaks:
		if elt[6] is False:
			# f1.write(str(elt)+"\n")
			f1.write(str(elt[0])+"\t"+str(elt[4])+"\t"+str(elt[5])+"\t"+str(elt[3])+"\n")

with open(data+file1+"_"+file2+"_peaks_merged.bed","w") as f1:
	for elt in all_peaks: # pour chaque peak
		if len(elt) > 7: # si merged
			# f1.write(str(elt)+"\n")
			mini=elt[1]
			maxi=elt[2]
			regions=elt[3]
			notcontained=False
			for elt2 in elt:
				if type(elt2) == list:
					regions="both" # str(regions)+","+str(elt2[3])
					if (int(elt2[6])==1): # elt2 contained
						f1.write(str(elt[0])+"\t"+str(elt[1])+"\t"+str(elt[2])+"\t"+str(regions)+"\n")
					else:
						notcontained=True
						if elt2[1]<mini:
							mini=elt2[1]
							if elt2[2]>maxi:
								maxi=elt2[2]
			if notcontained:
				f1.write(str(elt[0])+"\t"+str(mini)+"\t"+str(maxi)+"\t"+str(regions)+"\n")
		
