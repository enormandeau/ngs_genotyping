Get problematic individuals
- put problematic alleles in problematic_alleles.txt
- sort problematic_alleles.txt -V | uniq -c | sort -nr

Have a look at problematic alleles
c; grep -E "A15[[:space:]]|A15$" genotypes_jinquan_salmo_salar.tsv_output_table.csv | wc -l

Rince Later Repeat

