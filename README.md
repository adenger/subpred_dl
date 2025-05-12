# subpred_deeplearning

## create environment

mamba env create -f environment.yml
conda activate subpred_deeplearning

pip install -e .

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

Uniref version 2022_0? 

All except uniref updated on 11.05.2025

<!-- TODO backup and restore of pickles -->

<!-- TODO docker container with only data/datasets. -->

<!-- https://github.com/agemagician/ProtTrans/blob/master/Embedding/prott5_embedder.py -->
