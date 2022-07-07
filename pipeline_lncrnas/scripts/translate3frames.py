#!/usr/bin/env python3
"""
Script to translate sequences in the 3 possible frames
Arguments: 1/ imput:dna_sequences 2/ output:protein_sequences
"""
import sys
from Bio import SeqIO

translatedfile = open(sys.argv[2], "w")
for record in SeqIO.parse(sys.argv[1], "fasta"):
    t1 = record.seq.translate()
    translatedfile.write(">" + record.id + "\n")
    translatedfile.write(str(t1) + "\n")
    t2 = record[1:].seq.translate()
    translatedfile.write(">" + record.id + "\n")
    translatedfile.write(str(t2) + "\n")
    t3 = record[2:].seq.translate()
    translatedfile.write(">" + record.id + "\n")
    translatedfile.write(str(t3) + "\n")
    
translatedfile.close()
