#!/bin/bash

# After running mapreadstogenome_forme_HISAT2StringTie.sh, this script pulls TPM values from the abundance merged tab files and saves them in a CSV file. 

mkdir -p ./results/merged_TPM_all

while read name; do

	cp ./results/stringtie/${name}/${name}_abund_merged.tab ./results/merged_TPM_all

done < ./temp/basenamereads.txt

# Combine TPM values
# Columns 1=geneID, 2=gene_name, 3=reference, 4=strand, 5=start, 6=end, 7=coverage, 8=FPKM, 9=TPM
# I am collecting TPM, but you can collect as many columns as needed. Adjust -f parameter in line 21, and modify or replicate line 22 for each column.

while read name; do 

	echo $name
	cat ./results/merged_TPM_all/${name}_abund_merged.tab | tr "\\t" "," > ./results/merged_TPM_all/${name}.csv
	cut -d, -f 1,9 ./results/merged_TPM_all/${name}.csv > ./results/merged_TPM_all/${name}_values.csv
	sed -i "1s/TPM/$name\_TPM/" ./results/merged_TPM_all/${name}_values.csv

done < ./temp/basenamereads.txt 

paste -d , ./results/merged_TPM_all/*_values.csv > ./results/merged_TPM_all/all_values.csv

echo "Make sure you correctly filter your data for TPM of 0 values!"

