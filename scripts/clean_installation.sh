#!/bin/bash
# Remove files in working folders to create a fresh installation of ngs_genotyping
# WARNING: This will delete any work done in this ngs_genotyping installation!

echo "WARNING: this command will delete any work done in this ngs_genotyping installation!"
read -p "Are you sure you want do proceed? " -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm blast_database/* 2> /dev/null
    rm blast_results/* 2> /dev/null
    rm cleaned_aligned/* 2> /dev/null
    rm fasta/* 2> /dev/null
    rm individual_summary/* 2> /dev/null
    rm raw_data/* 2> /dev/null
    rm trimmed_separated_sequences/* 2> /dev/null
    echo -e "\n--All working folders have been emptied--"
fi

