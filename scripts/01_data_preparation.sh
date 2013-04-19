#!/bin/bash

MIDNAMES="trimmed_separated_sequences/midnames.txt"

# Cleanup
rm raw_data/*identified* 2> /dev/null

# Renaming the sequences with the MIDs
echo "Renaming sequences..."
echo "" > $MIDNAMES  # Create file containing the names of the MIDs
for file in `ls -1 raw_data/*MID*.fna | sed 's/fna//'`; do  \
    mid=`echo $file | perl -pe 's/.*([0-9]{2})\.(MID[0-9]*).*/\2-\1/'`; \
    echo $mid >> $MIDNAMES
    perl -pe 's/>/>@ARGV[0]_/;' $file"fna" $mid >$file"fna.identified" 2> /dev/null; \
    perl -pe 's/>/>@ARGV[0]_/' $file"qual" $mid >$file"qual.identified" 2> /dev/null; \
done
perl -i -pe 's/^\n//' $MIDNAMES

# Concatenate fasta and qual files after renaming
cat raw_data/*.fna.identified > raw_data/all_identified.fna
cat raw_data/*.qual.identified > raw_data/all_identified.qual

#gives the frequency distribution of sequence lenghts for a given fasta file.
./scripts/fasta_unwrap.py raw_data/all_identified.fna raw_data/all_identified.fna.unwrap
grep -v '>' raw_data/all_identified.fna.unwrap | \
    awk '{print length()}' | \
    sort -n | \
    uniq -c | \
    awk '{ print $2, $1}' > raw_data/all_identified_length_hist.txt
rm raw_data/all_identified.fna.unwrap

# Plot graph of sequence lengths
echo "
Examine the distribution of sequence lengths from the gnuplot graph to decide on
the minimal and maximal length thresholds that should be applied to remove reads
that are either too short or too long.

If you want to look at the graph again without running the whole script, type
the followint command in the terminal:

    ./scripts/plotlines.sh ./raw_data/all_identified_length_hist.txt
"
./scripts/plotlines.sh ./raw_data/all_identified_length_hist.txt

