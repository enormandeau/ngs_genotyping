#!/bin/bash

# Global variables
TRIM_LEFT=8 # The sequences should already be trimmed appropriately
TRIM_RIGHT=35 # idem
INDELS=$1
SNPS=$2
TEMP_FILE="clean_align_clean.temp"

# Help
echo "----------------------------------------------------------"
echo "Help for program 05_align_clean_align.sh"
echo
echo "Usage:"
echo "  ./scripts/05_align_clean_align.sh INDELS SNPS"
echo
echo "  INDELS = Decimal number, minimal proportion for indels to be considered real"
echo "  SNPS = Decimal number, minimal proportion for SNPs to be considered real"
echo
echo "Example run:"
echo "  ./scripts/05_align_clean_align.sh 0.01 0.01"
echo "----------------------------------------------------------"

# Create list of files to treat
ls -1 trimmed_separated_sequences/*MID*.fasta > $TEMP_FILE
sed -i 's/\.fasta//' $TEMP_FILE

# Clean, align, Clean, align, Clean
cat $TEMP_FILE | while read i; do
    echo "Treating file: $i"
    muscle -in "$i.fasta" -out $i".aca_temp1" -diags -maxiters 4 -quiet
    ./scripts/alignment_clean.py -i $i".aca_temp1" -o $i".aca_temp2" -l $TRIM_LEFT -r $TRIM_RIGHT -I $INDELS -s $SNPS
    muscle -in $i".aca_temp2" -out $i".aca_temp3" -diags -maxiters 4 -quiet
    ./scripts/alignment_clean.py -i $i".aca_temp3" -o $i".aca_temp4" -I 0.025 -s 0.025
    muscle -in $i".aca_temp4" -out $i"_cleaned_aligned.fasta" -diags -maxiters 4 -quiet
    rm trimmed_separated_sequences/*aca_temp*
done

# Clean up directory
rm $TEMP_FILE
mv trimmed_separated_sequences/*_cleaned_aligned.fasta cleaned_aligned

