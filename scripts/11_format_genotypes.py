#!/usr/bin/python
"""
Format the genotypes obtained from the MHC pipeline to a usefull format

Usage:
    program.py <genotypes.txt> <annotation.txt>

The <genotypes.txt> file contains the output of phase two from Higgy Pop

The <annotation.txt> file is an optional tab separated file containing one MID
per line in the first column, followed by any pertinent annotation data to be
joined with the genotypes. It is essential that the MID format be exactly the
same as found in the 'individual_summaries_combined' file. Eg: MID125-04.

"""

# Importing modules
import sys
from collections import defaultdict

# Parsing command line
try:
    genotype_file = sys.argv[1]
except:
    print __doc__
    sys.exit(1)

try:
    annotation_file = sys.argv[2]
except:
    annotation_file = None

genotype_stub = genotype_file.replace(".txt", "")

# Loading file
data_ind = [l.strip().split("_") for l in open(genotype_file).readlines() if l.strip() != ""]
try:
    data_annot = [l.strip().split("\t") for l in open(annotation_file).readlines() if l.strip() != ""]
except:
    data_annot = None

# Creating dictionary of alleles in the individuals
individuals = defaultdict(set)
for i in data_ind:
    individuals[i[0]].add(i[2])

# Loading the individual annotation file
annotation = {}
if data_annot != None:
    for a in data_annot:
        annotation[a[0]] = ["\t".join(a)]
else:
    for i in data_ind:
        annotation[i[0]] = ["no annotation"]

# Generating first output: One individual-allele combination per line
with open(genotype_stub + "_output_list.csv", "w") as f:
    for i in sorted(individuals):
        for allele in sorted(individuals[i]):
            try:
                f.write(i + "\t" + annotation[i][0] + "\t" + str(len(individuals[i])) + "\t" + allele + "\n")
            except:
                f.write(i + "\t" + str(len(individuals[i])) + "\t" + allele + "\n")

# Generating second output: Table of individuals and alleles
with open(genotype_stub + "_output_table.csv", "w") as f:
    for i in sorted(individuals):
        try:
            f.write(i + "\t" + str(len(individuals[i])) + "\t" + \
                annotation[i][0] + "\t" + "\t".join(sorted(individuals[i])) + "\n")
        except:
            f.write(i + "\t" + str(len(individuals[i])) + "\t" + \
                "\t" + "\t".join(sorted(individuals[i])) + "\n")

