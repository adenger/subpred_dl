# subpred_deeplearning

## Installation

### Requirements

Linux, miniconda

### create environment

mamba env create -f environment.yml
conda activate subpred_deeplearning

### install code into environment

pip install -e .

## Recreating results from manuscript

Links to the pre-processed datasets used in the study are provided here:

data/datasets (..GB)

To re-create the results from pre-processed data, only the latter archive is needed. 

## Recreate pre-processing

The raw data can be downloaded here:

data/raw (..GB)

Uniref is version 2022_01 (contains enough proteins to create evolutionary profiles, and we already had pre-calculated PSSMs for most proteins from a previous project)

TODO place Uniref2022 in raw data folder.

The preprocessing of raw data can be carried out manually via

./preprocessing/create_blastdbs.sh

./preprocessing/create_datasets.py

Alternatively, a link to pre-processed data is provided above.

## Download/update and process raw data

To re-create the project entirely from scratch, the download and the initial pre-processing is also listed in the preprocessing/download*.sh scripts.

File formats of GO and Uniprot can change in the future, making the code incompatible. Changes between current and future database versions might affect the results.

Link to raw data is provided above.

All except Uniref updated on 11.05.2025

### download go annotations

./preprocessing/download_goa.sh

GOA UniProt (version 226), released on 06 May, 2025 and assembled using the publicly released data available in the source databases on 28 April, 2025.

### download uniprot data (conda env must be active)

./preprocessing/download_uniprot.sh

Uniprot Version 2025_02 released on 23.04.2025

### download go ontology

./preprocessing/download_go.sh

GO version 2025-03-16

### download uniref, create blast databases (conda env must be active)

./preprocessing/download_uniref.sh

Uniref version 2022_1

<!-- TODO docker container with only data/datasets. -->

<!-- https://github.com/agemagician/ProtTrans/blob/master/Embedding/prott5_embedder.py -->
