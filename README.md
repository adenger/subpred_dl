# subpred_deeplearning

## Installation

### Requirements

- Linux (tested on Ubuntu 24.04 LTS)
- [miniforge](https://github.com/conda-forge/miniforge)
- At least 64GB of memory
- Up to 500GB disk storage
- To re-calculate embeddings (optional): GPU compatible with CUDA Toolkit 12.6+

### Create environment, install code

```bash
conda env create -f environment_full.yml
conda activate subpred_deeplearning
pip install -e .
```

A separate environment is used for the DNN notebooks. All package versions are the same, except that the dnn_cpu environment uses the CPU version of tensorflow to train DNNs, whereas the subpred_deeplearning environment contains the CUDA-accellerated variant of the package to generate embeddings on the GPU. The DNN training is currently incompatible with the latest generation of Nvidia GPUs. The CPU version of TF also ensures full reproducibility of all results.

The SVM notebooks need to be executed first to generate the feature data, then the environment for the dnn notebooks (and the *plots_tables_human* notebook) can be recreated with:

```bash
conda env create -f environment_dnn_cpu_full.yml
conda activate dnn_cpu
pip install -e .
```

## Recreating results from manuscript

Links to the pre-processed datasets used in the study are provided here:

data/datasets (..GB)

To re-create the manuscript results from pre-processed data, only this archive is needed.

## Recreate pre-processing

The raw data can be downloaded here:

data/raw (..GB)

Uniref is version 2022_01 (contains enough proteins to create evolutionary profiles, and we already had pre-calculated PSSMs for most proteins from a previous project), everything else was downloaded on 11.05.2025.

TODO place Uniref2022 in raw data folder.

The preprocessing of raw data can then be carried out manually via:

TODO turn notebook into py

```bash
./preprocessing/create_blastdbs.sh
./preprocessing/create_3Di_fasta.sh
./preprocessing/create_datasets.py  # TODO
```

## Download/update and process raw data

To re-create the project entirely from scratch, the download and the initial pre-processing is also listed in the preprocessing/download*.sh scripts.

File formats of GO and Uniprot can change in the future, making the code incompatible. Changes between current and future database versions might affect the results.

Link to raw data is provided above.

All except Uniref updated on 11.05.2025

### Download go annotations

```bash
./preprocessing/download_goa.sh
```

GOA UniProt (version 226), released on 06 May, 2025 and assembled using the publicly released data available in the source databases on 28 April, 2025.

### Download AlphafoldDB PDB files for model organisms

```bash
./preprocessing/download_alphafolddb.sh
```

Version 4 from 2022, downloaded on 15.05.2025

Additional tar files from alphafolddb (https://www.alphafold.ebi.ac.uk/download) can be added to the script, to include more organisms. They will automatically be pre-processed.

### Download uniprot data (conda env must be active)

```bash
./preprocessing/download_uniprot.sh
```

Uniprot Version 2025_02 released on 23.04.2025

### Download go ontology

```bash
./preprocessing/download_go.sh
```

GO version 2025-03-16

### Download uniref, create blast databases (conda env must be active)

```bash
./preprocessing/download_uniref.sh
```

Uniref version 2022_1

### Download interpro annotation names

```bash
./preprocessing/download_interpro.sh
```

Downloaded version 2025-04-22

<!-- TODO docker container with only data/datasets. -->

<!-- https://github.com/agemagician/ProtTrans/blob/master/Embedding/prott5_embedder.py -->
