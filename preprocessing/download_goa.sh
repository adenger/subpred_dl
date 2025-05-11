#!/bin/bash
echo "downloading and filtering uniprot go annotations..."
wget https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz -O - | gunzip -c | awk 'BEGIN {OFS="\t";FS="\t"} ($1 == "UniProtKB") {print $2,$4,$5,$7,$9,$14}' | sort -u | gzip > data/raw/goa/goa_uniprot_all_ebi_filtered.tsv.gz
echo "done"
