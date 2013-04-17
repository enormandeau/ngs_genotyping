#!/usr/bin/python
"""Use files with frequency counts of blasts on possible alleles (see format
below) to genotype individuals.

Usage:
    genotype_from_blast_results.py  file_list  threshold_file  output_file

- file_list is a file containing the path to the individual summary files. the
file_list format is as follows:

individual_summary/MID100-02_cleaned_aligned.fasta.blasts_summary
individual_summary/MID100-03_cleaned_aligned.fasta.blasts_summary
individual_summary/MID100-04_cleaned_aligned.fasta.blasts_summary
...
individual_summary/MID102-03_cleaned_aligned.fasta.blasts_summary

- The individual summary files are in the following format (output format of
sort | uniq -c | sort -nr in Linux terminal):

    316 A_27
    197 A_15
     13 A_1
      4 A_8
      4 A_19
      3 A_28
      1 A_9
      1 A_16

- threshold_file is a file with, on each line, the name of an allele and the
minimal count to consider this allele real in an individual, separated by a
tabulation. Format:

0.75 A_1
0.23 A_2
0.10 A_7
0.10 A_10

- output_file is the name of the file that will contain the results
"""

# Importing modules
import sys
import re

# Main
if __name__ == '__main__':
    try:
        file_list = sys.argv[1]
        threshold_file = sys.argv[2]
        output_file = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)
    
    to_write = []
    individuals = [l.strip() for l in open(file_list) if l.strip() != ""]
    allele_threshold ={}
    thresholds = [l.strip().split("\t") for l in open(threshold_file) if l.strip() != ""]
    
    for tr in thresholds:
        t, a = float(tr[0]), tr[1]
        allele_threshold[a] = t
    
    for i in individuals:
        itemp = i.replace("individual_summary/", "")
        itemp = itemp.replace("A_", "")
        itemp = re.sub('(MID[0-9]+-[0-9]+).*', '\\1', itemp)
        data = [l.strip() for l in open(i) if l.strip() != ""]
        allele_dict = {}
        for d in data:
            count, allele = d.split(" ")
            allele_dict[allele] = float(count)
        total = sum(allele_dict.values())
        allele_number = 0
        #temp_min_proportion = min_proportion
        for a in allele_dict:
            atemp = a.replace("_", "")
            if allele_dict[a] >= allele_threshold[a]:
                allele_number += 1
                to_write.append("_".join([itemp, str(allele_number), atemp]))
    open(output_file, "w").write("\n".join(to_write))



