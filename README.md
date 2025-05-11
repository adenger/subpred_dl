# subpred_deeplearning

## create environment

mamba env create -f environment.yml
conda activate subpred_deeplearning

## download go annotations

./preprocessing/download_goa.sh

GOA UniProt (version 226), released on 06 May, 2025 and assembled using the publicly released data available in the source databases on 28 April, 2025.

## download uniprot data (conda env must be active)

./preprocessing/download_uniprot.sh

Uniprot Version 2025_02 released on 23.04.2025

## download go ontology

./preprocessing/download_go.sh

GO version 2025-03-16

## download uniref, create blast databases (conda env must be active)

./preprocessing/download_uniref.sh

./preprocessing/create_blastdbs.sh

Uniref version 2025_02 released on 23.04.2025

All downloaded on 11.05.2025