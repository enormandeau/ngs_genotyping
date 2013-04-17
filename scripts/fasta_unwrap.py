#!/usr/bin/python
import sys
from Bio import SeqIO

in_file = open(sys.argv[1], "rU")
out_file = open(sys.argv[2], "w")
sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(in_file, 'fasta'))
with open(sys.argv[2], "w") as out_file:
    for seq in sequences:
        out_file.write(">" + seq[0] + "\n" + seq[1] + "\n")
