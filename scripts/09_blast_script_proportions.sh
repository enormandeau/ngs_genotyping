#!/bin/bash
# Create blast database
evalue=$1
maxnum=`grep ">" allele_database.fasta | sed -E 's/[^0-9]*//g' | sort -n | tail -1`
database_file="allele_database.fasta"
database_title="allele_database"
sed -i -E 's/\-//g' $database_file

# Empty result folders
rm ./blast_results/* 2> /dev/null
rm ./individual_summary/* 2> /dev/null

# Remove dashes '-' from fasta sequences
sed -i -E 's/\-//g' *.fasta

# Create blast database
makeblastdb -in $database_file -title $database_title -dbtype nucl -out ./blast_database/$database_title

# Blast sequences
for file in `ls -1 fasta/*.fasta | sed -E 's/^.+\///'`; do blastn -db ./blast_database/$database_title -query fasta/$file -out blast_results/$file".blast" -max_target_seqs 1 -num_threads 8 -outfmt 6 -evalue 1e-$evalue; awk  '{print $2}' blast_results/$file".blast" | sort | uniq -c | sort -nr >individual_summary/$file".blasts_summary"; done

# Create blast summaries
for file in `ls -1 blast_results/* | sed -E 's/^.+\///'`; do awk 'BEGIN {FS="\t"} {print $2}' blast_results/$file; done | sort | uniq -c | sort -nr > summary_all.txt

# Create threshold file
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read i; do echo -e "0.1\t$i"; done > individual_thresholds.txt

# Figures of the distribution of evalues of the blasts on each allele
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read i; do cat ./blast_results/*.blast | grep "$i[^0-9]" | awk '{print $11}' | sed -E 's/[0-9]e\-//' | sort -g | uniq -c > ./blast_results/data_$i; ./scripts/plot2lines.sh ./blast_results/data_$i; done

# Calculate proportions
for file in $(ls -1 ./individual_summary/*summary); do echo Treating file $file; tot=$(perl -pe 's/^ +//' $file | cut -d " " -f 1 | ./scripts/total.sh); cat $file | while read line; do value=$(echo $line | perl -pe 's/^ +//' | cut -d " " -f 1); allele=$(echo $line | perl -pe 's/^ +//' | cut -d " " -f 2); prop=$(echo "scale=9; $value / $tot" | bc | cut -c -7); echo $prop $allele; done > $file".proportions"; done

# Figures of the number of sequences per allele in individuals where it is found
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read i; do grep "$i$" ./individual_summary/*summary.proportions | perl -pe 's/^.+://' | sort -nr | awk '{print $1}' > ./individual_summary/data_$i; ./scripts/plot.sh individual_summary/data_$i; done

# Explore depth of sequences per individual
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read allele; do echo START $allele; for ind in $(grep " $allele$" ./individual_summary/*summary.proportions | perl -pe 's/:/ /' | sort -t " " -k 2 -nr | perl -pe 's/ .+//'); do echo $ind; cat $ind; done; done > individual_summaries_combined.txt

# Create list of filenames for the .proportions files
ls -1 individual_summary/*summary.proportions > summary_file_names_proportions.txt

