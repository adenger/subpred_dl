#!/bin/bash
wget -P data/raw/alphafolddb/tarfiles https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000006548_3702_ARATH_v4.tar
wget -P data/raw/alphafolddb/tarfiles https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000625_83333_ECOLI_v4.tar
wget -P data/raw/alphafolddb/tarfiles https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
wget -P data/raw/alphafolddb/tarfiles https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar

# additional tar files from alphafolddb can be added here, to add more proteins (https://www.alphafold.ebi.ac.uk/download)