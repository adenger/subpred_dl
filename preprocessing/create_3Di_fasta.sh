#!/bin/bash
for file in data/raw/alphafolddb/tarfiles/*.tar; do tar xf "$file" --directory=data/raw/alphafolddb/pdbs --wildcards "*.pdb.gz"; done
foldseek createdb data/raw/alphafolddb/pdbs data/raw/alphafolddb/foldseekdb/queryDB
# see https://github.com/steineggerlab/foldseek/issues/15#issuecomment-1065876787
foldseek lndb data/raw/alphafolddb/foldseekdb/queryDB_h data/raw/alphafolddb/foldseekdb/queryDB_ss_h
foldseek convert2fasta data/raw/alphafolddb/foldseekdb/queryDB_ss data/raw/alphafolddb/3Di_sequences.fasta
