# -*- coding: utf-8 -*-

print("")
print("====> Launching merge_all_peaksV2.py from bin/" )
print("")

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
		l[1]=int(l[1]) # start
		l[2]=int(l[2]) # end
		l.append(False) # = l[4]
		all_peaks.append(l) # l -> chr start end name False
   

	 
i=0
j=0
k=0
while (i < len(all_peaks) -1):
	j=i+1
	
	while (j < len(all_peaks) -1): # dual loop for sliding on the lines / start J >= start I
		# Start J >= start I
		if all_peaks[i][0] == all_peaks[j][0] and all_peaks[i][4] is not True : # même CHR # TODO enlever condition True ?
			if  all_peaks[j][2] <= all_peaks[i][2] : # stop de J avant stop I / I contient J
				all_peaks[j][4]=True # a.k.a. merged
				all_peaks[i][4]=True # a.k.a. merged
				all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3]]) 
				if (all_peaks[j][1] - all_peaks[i][1]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 : 
					# If length of I before start of J higher than 50% of the peak -> create a new spe peaks out of this - not re tested for merge !
					all_peaks.insert(i,[all_peaks[i][0],all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],False])
					i+=1
					j+=1
				if (all_peaks[i][2] - all_peaks[j][2]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 :
					# If length of I after end of J higher than 50% of the peak -> create a new spe peaks out of this  - re tested for potential merge
					all_peaks.insert(j+1,[all_peaks[i][0],all_peaks[j][2]+1,all_peaks[i][2],all_peaks[i][3],False])
			elif ((all_peaks[j][1] < all_peaks[i][2]) and (all_peaks[j][2] > all_peaks[i][2])) : # start de J contenu dans I mais pas le stop.
				if all_peaks[i][2] -all_peaks[i][1] < all_peaks[j][2] -all_peaks[j][1]: # I plus petit que J
					if (all_peaks[i][2] - all_peaks[j][1]) /float(all_peaks[i][2] - all_peaks[i][1]) > 0.8:
						all_peaks[j][4]=True # a.k.a. merged
						all_peaks[i][4]=True # a.k.a. merged
						all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3]]) 
						if (all_peaks[j][2] - all_peaks[i][2]) / float(all_peaks[j][2] - all_peaks[j][1] ) > 0.5 : # Non overlap de J / taille de J > 0.5 
							all_peaks.insert(j+1,[all_peaks[i][0],all_peaks[i][2]+1,all_peaks[j][2],all_peaks[j][3],False])
				else: # I plus grand ou égal à J
					if (all_peaks[i][2] - all_peaks[j][1]) /float(all_peaks[j][2] - all_peaks[j][1]) > 0.8: # overlap / taille de J > 0.8
						all_peaks[j][4]=True # a.k.a. merged
						all_peaks[i][4]=True # a.k.a. merged
						all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3]])
						if (all_peaks[j][1] - all_peaks[i][1]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 : # Non overlap de I / taille de I > 0.5
							all_peaks.insert(i,[all_peaks[i][0],all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],False])
							i+=1
							j+=1						
			elif ((all_peaks[j][1] >= all_peaks[i][2])): #JS après IE -> pas de merge, on peut arreter le while il n'y aura plus de peak à merge après
				k+=1
				break
			# la
		else : # CHR différents, plus de merge possible après
			break
		# if j-i == 2:
		#	 print j-i,all_peaks[i][0],all_peaks[i][1],all_peaks[i][2],all_peaks[j][1],all_peaks[j][2]
		j+=1
	i+=1
	if  i%1000 ==0 :
		 print(i)
		
print("k="+str(k))
print("i="+str(i))
print("j="+str(j))


with open(output+"/"+file1+"_"+file2+"_peaks_uniques.bed","w") as f1:
	for elt in all_peaks:
		if elt[4] is False:
			f1.write(str(elt[0])+"\t"+str(elt[1])+"\t"+str(elt[2])+"\t"+str(elt[3])+"\n")

with open(output+"/"+file1+"_"+file2+"_peaks_merged.bed","w") as f1:
	for elt in all_peaks: # pour chaque peak
		if len(elt) > 5: # si merged
			mini=elt[1]
			maxi=elt[2]
			regions=elt[3]
			
			for elt2 in elt:
				if type(elt2) == list:
					regions=[str(regions),str(elt2[3])]
					regions.sort()
					regions="&".join(regions)
					#regions=str(regions)+"&"+str(elt2[3]) #"both"
					if elt2[1]>mini:
						mini=elt2[1]
						if elt2[2]<maxi:
							maxi=elt2[2]
			f1.write(str(elt[0])+"\t"+str(mini)+"\t"+str(maxi)+"\t"+str(regions)+"\n")
		

