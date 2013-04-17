#!/bin/bash
# Create blast database
evalue="126"
maxnum=`grep ">" allele_database.fasta | sed -E 's/[^0-9]*//g' | sort -n | tail -1`
database_file="allele_database.fasta"
database_title="allele_database"
sed -i -E 's/\-//g' $database_file

# Empty result folders
rm ./blast_results/*
rm ./individual_summary/*

# Remove dashes '-' from fasta sequences
sed -i -E 's/\-//g' *.fasta

# Create blast database
makeblastdb -in $database_file -title $database_title -dbtype nucl -out ./blast_database/$database_title

# Blast sequences
for file in `ls -1 fasta/*.fasta | sed -E 's/^.+\///'`; do blastn -db ./blast_database/$database_title -query fasta/$file -out blast_results/$file".blast" -max_target_seqs 1 -num_threads 8 -outfmt 6 -evalue 1e-$evalue; awk  '{print $2}' blast_results/$file".blast" | sort | uniq -c | sort -nr >individual_summary/$file".blasts_summary"; done

# Create blast summaries
for file in `ls -1 blast_results/* | sed -E 's/^.+\///'`; do awk 'BEGIN {FS="\t"} {print $2}' blast_results/$file; done | sort | uniq -c | sort -nr > summary_all.txt

# Create threshold file
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read i; do echo -e "2\t$i"; done > individual_thresholds.txt

# Figures of the distribution of evalues of the blasts on each allele
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read i; do cat ./blast_results/*.blast | grep "$i[^0-9]" | awk '{print $11}' | sed -E 's/[0-9]e\-//' | sort -g | uniq -c >./blast_results/data_$i; ./scripts/plot2lines.sh ./blast_results/data_$i; done

# Figures of the number of sequences per allele in individuals where it is found
grep ">" allele_database.fasta | sed -E 's/>(A_[0-9]+).+?$/\1/; s/>//' | while read i; do grep "$i$" ./individual_summary/*summary | sort -t " " -k 2 -nr | awk '{print $2}' > ./individual_summary/data_$i; ./scripts/plot.sh ./individual_summary/data_$i; done

# Explore depth of sequences per individual
for i in `grep "A_8$" ./individual_summary/*summary | sort -t " " -k 2 -nr | awk '{print $1}' | sed -E 's/://'`; do echo $i; cat $i; done > individual_summaries_combined.txt

grep -E ' A_48$' individual_summaries_combined.txt

# Create file listing the summary files
ls -1 individual_summary/*summary > summary_file_names.txt

