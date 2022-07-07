#!/usr/bin/env python3
"""
Script to filter the gffcompare output for novel potential lncRNAs
arguments= 1) gtf file from gffcompare 2) output directory(/some/directory)
created output files:
iux.gtf(u,i,x) = table of i, u and x class transcripts
potential_lncrna.gtf = gtf with filtered transcripts
potentialids.txt = list of ids of filtered transcripts
"""
import re
import sys

#separate transcripts by class and get ids in list
other = []
with open(sys.argv[1]) as f:
    for line in f:
        if re.search("class_code \"i\"", line):
            tabline = re.split("\t|; ", line)
            tabline[8] = tabline[8][15:len(tabline[8])-1]
            other.append(tabline[8])
        elif re.search("class_code \"u\"", line):
            tabline = re.split("\t|; ", line)
            tabline[8] = tabline[8][15:len(tabline[8])-1]
            other.append(tabline[8])
        elif re.search("class_code \"x\"", line):
            tabline = re.split("\t|; ", line)
            tabline[8] = tabline[8][15:len(tabline[8])-1]
            other.append(tabline[8])

#get lncRNA transcripts and exons (i, u and x classes) from gtf file
iuxgtf = open(sys.argv[2] + "/iux.gtf", "w")

with open(sys.argv[1]) as f:
    for line in f:
        tabline = re.split("\t|; ", line.rstrip())
        myid = tabline[8][15:len(tabline[8])-1]
        tabline[8] = tabline[8][15:len(tabline[8])-1]
        tabline[9] = tabline[9][9:len(tabline[9])-1]
        if myid in other:
            iuxgtf.write("\t".join(tabline) + "\n")

iuxgtf.close()
print("Classes i, u and x stored in iux.gtf")

#remove single exon transcripts < 500 and multi exon transcripts < 200
iuxgtf_list = []
with open(sys.argv[2] + "/iux.gtf") as f:
    for line in f:
        tabline = re.split("\t", line.rstrip())
        lenght = int(tabline[4]) - int(tabline[3])
        tabline.append(lenght)
        iuxgtf_list.append(tabline)

singleexon = []
multiexon = []
for entry in iuxgtf_list:
    if entry[2] == "transcript":
        exons = []
        for entry2 in iuxgtf_list:
            if entry2[8] == entry[8] and entry2[2] == "exon":
                exons.append(entry2[-1])
        if (len(exons) == 1) and (10000 > exons[0] > 500):
            singleexon.append(entry)
        elif (len(exons) > 1) and (50000 > sum(exons) > 200):
            multiexon.append(entry)

print("Removed single exon transcripts < 500, > 10000 and multi exon transcripts < 200")

#Save single and multi exon transcripts in correct gtf file
potentialidstxt = open(sys.argv[2]+ "/potentialids.txt", "w")
potentialids = []
for entry in singleexon:
    potentialids.append(entry[8])
for entry in multiexon:
    potentialids.append(entry[8])
for entry in potentialids:
    potentialidstxt.write(entry + "\n")
potentialidstxt.close()
print("Potential lnc ids saved in potentialids.txt")

print("Searching gtf file for potential lncrna ids:")
potentialgtf = open(sys.argv[2] + "/potential_lncrna.gtf", "w")

with open(sys.argv[1]) as f:
    for entry in f:
        tabline = re.split("\t|; ", entry)
        myid = tabline[8][15:len(tabline[8])-1]
        if myid in potentialids:
            potentialgtf.write(entry)

potentialgtf.close()
print("Finished: potential lncRNAs added to gtf file")
