#!/bin/bash
wget -P data/raw/alphafolddb https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000006548_3702_ARATH_v4.tar
wget -P data/raw/alphafolddb https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000625_83333_ECOLI_v4.tar
wget -P data/raw/alphafolddb https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
wget -P data/raw/alphafolddb https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar

tar xf data/raw/alphafolddb/UP000006548_3702_ARATH_v4.tar --directory=data/datasets/pdb/ --wildcards "*.pdb.gz"
tar xf data/raw/alphafolddb/UP000000625_83333_ECOLI_v4.tar --directory=data/datasets/pdb/ --wildcards "*.pdb.gz"
tar xf data/raw/alphafolddb/UP000005640_9606_HUMAN_v4.tar --directory=data/datasets/pdb/ --wildcards "*.pdb.gz"
tar xf data/raw/alphafolddb/UP000002311_559292_YEAST_v4.tar --directory=data/datasets/pdb/ --wildcards "*.pdb.gz"