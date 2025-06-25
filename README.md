# Membrane transporter substrate prediction with deep learning methods

## Installation

### Requirements

- Linux (tested on Ubuntu 24.04 LTS in WSL2)
- [miniforge](https://github.com/conda-forge/miniforge)
- To recreate feature datasets (optional): 
    - Up to 500GB disk storage (for BLAST databases), 64GB Memory (psiblast)
    - GPU compatible with CUDA Toolkit 12.6+, >=16GB VRAM (embeddings)

### Create environment, install project code as python package

```bash
conda env create -f environment_full.yml
conda activate subpred_deeplearning
pip install -e .
```

A separate environment is used for the DNN notebooks (notebooks 08-14). All package versions are the same, except that the dnn_cpu environment uses the CPU version of tensorflow to train DNNs, whereas the subpred_deeplearning environment contains the CUDA-accellerated variant of the package to generate embeddings on the GPU. The DNN training is currently incompatible with the latest generation of Nvidia GPUs. The CPU version of TF also ensures full reproducibility of all results.

```bash
conda env create -f environment_dnn_cpu_full.yml
conda activate dnn_cpu
pip install -e .
```

## Recreating results from manuscript

The raw data is available here:

[/data/raw (113GB)](https://1drv.ms/u/c/886666fa46e5db95/EYpQvp0O-G9EsvFw4ZFaSvwB3LSWCF0L2qD77VUccFlUnQ?e=82mHc5)

Running the **01_preprocessing** notebook will turn the raw data into pre-processed pickles. To vastly speed up the feature computation, we saved the PSSMs and embeddings that we calculated for all proteins in the dataset in a cache folder. Once they are extracted into the appropriate folder, the feature generation methods will read these files instead of calculating everything from scratch. The preprocessed pickles, along with cached PSSMs and embeddings, are available for download here:

[/data/datasets (1.4GB)](https://1drv.ms/u/c/886666fa46e5db95/EaTz162K0i9AlGv-kczsR44Bgc4pXxC2i4OyMkJ0VhwUuQ?e=F9UKgM)

After extracting the data into the matching folders (tar -xf from the root directory of the repository), the notebooks can be re-calculated. Here, it is important to run the svm notebooks (02-07) first with the *subpred_deeplearning* conda environment, and then the dnn notebooks (08-14) with the *dnn_cpu* environment, for reasons mentioned above. The ML feature data that is created by the SVM notebooks and subsequently read by the DNN notebooks is, alternatively, also saved in an archive that can be downloaded here:

[/data/tmp_data (37MB)](https://1drv.ms/u/c/886666fa46e5db95/EUXwykta-7pLsUj1dTIhfk0B8VptAZYJ-RHs3LyHNYWueg?e=EukKna)

Finally, the evaluation scores from all iterations of the repeated 5-fold cross validation, along with generated plots, are available here:

[/data/results (3MB)](https://1drv.ms/u/c/886666fa46e5db95/ESccd6GL03lLsKx_pHsF3KcBVH4qUwkp15f2E04ffg6vtA?e=RDzgEU)

## How the raw data was assembled (do not run)

All commands used to assemble /data/raw were saved in the preprocessing folder. **Note that these scripts always download the latest version of each database, and the contents of the datasets might change in the future**. Uniref is version 2022_01 (contains enough proteins to create evolutionary profiles, and we already had pre-calculated PSSMs for most proteins from a previous project), everything else was downloaded on 11.05.2025. They were executed in this order:

### GO annotations

```bash
./preprocessing/download_goa.sh
```

GOA UniProt (version 226), released on 06 May, 2025 and assembled using the publicly released data available in the source databases on 28 April, 2025.

### AlphafoldDB PDB files for model organisms

```bash
./preprocessing/download_alphafolddb.sh
```

Version 4 released in 2022, downloaded on 15.05.2025

Additional tar files from alphafolddb (https://www.alphafold.ebi.ac.uk/download) can be added to the script, to include more organisms. They will automatically be pre-processed.

### Uniprot data

```bash
conda activate subpred_deeplearning
./preprocessing/download_uniprot.sh
```

Uniprot Version 2025_02 released on 23.04.2025

### GO OBO

```bash
./preprocessing/download_go.sh
```

GO version 2025-03-16

### Uniref

```bash
./preprocessing/download_uniref.sh
```

Uniref version 2022_1

### Interpro annotation names

```bash
./preprocessing/download_interpro.sh
```

Downloaded version 2025-04-22

### Create 3Di fasta file and blast databases:

```bash
conda activate subpred_deeplearning
./preprocessing/create_blastdbs.sh
./preprocessing/create_3Di_fasta.sh
```


<!-- TODO docker container with only data/datasets. -->

<!-- https://github.com/agemagician/ProtTrans/blob/master/Embedding/prott5_embedder.py -->
