#!/usr/bin/env python
"""
Extract list of non-coding transcripts for each program and make venn diagram
Also filter out lncRNAs with more than 10 transcripts
Arguments: 
    1/ cpat results 
    2/ feelnc results
    3/ hmmer results
    4/ output: lnc id list
    5/ output: venn diagram
"""
# extract list of non-coding genes transcripts for each program and make venn diagram
import matplotlib
matplotlib.use('Agg')
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import re
import sys
import pandas as pd

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(sys.argv[5])

cpat_nc = []
feelnc_nc = []
hmmer_nc = []
all = []
hmmer_c = []

with open(sys.argv[1]) as f:
    next(f)
    for line in f:
        tabline = re.split("\t", line.rstrip())
        all.append(tabline[0])
        if len(tabline) > 1:
            if float(tabline[10]) < 0.4:
                cpat_nc.append(tabline[0])
        else:
            cpat_nc.append(tabline[0])

with open(sys.argv[2]) as f:
    for line in f:
        if re.search("transcript", line):
            tabline = re.split("\t|; ", line.rstrip())
            mstrg = tabline[8][15:len(tabline[8])-1]
            feelnc_nc.append(mstrg)

with open(sys.argv[3]) as f:
    for line in f:
    	tabline = line.rstrip().split()
    	if tabline[0] != "#":
    		if tabline[2] not in hmmer_c:
    			hmmer_c.append(tabline[2])

for entry in all:
	if entry not in hmmer_c:
		hmmer_nc.append(entry)

print("all " + str(len(all)))
print("cpat " + str(len(cpat_nc)))
print("feelnc " + str(len(feelnc_nc)))
print("hmmer " + str(len(hmmer_nc)))

venn3([set(cpat_nc), set(feelnc_nc), set(hmmer_nc)], set_labels = ('CPAT', 'FEELnc', 'HMMER'))
plt.savefig(sys.argv[5])

# Remove lncRNA genes with many isoforms (>10) and write lnc id list file
allunfilter = []
for entry in all:
    if entry in cpat_nc and entry in feelnc_nc and entry in hmmer_nc:
        tabline = entry.rstrip().split(".")
        tabline.append(entry.rstrip())
        allunfilter.append(tabline)

allunfilter_df = pd.DataFrame(allunfilter)

allfilter = []
for entry in allunfilter:
    if len(allunfilter_df[allunfilter_df[1] == entry[1]]) < 11:
        allfilter.append(entry[3])

lncrna = open(sys.argv[4], "w")
for entry in allfilter:
	lncrna.write(entry + "\n")
lncrna.close()
