#!/bin/bash
echo "Extracting Uniref50..."
gunzip -c data/raw/uniref/uniref50/uniref50.fasta.gz > data/datasets/blastdb/uniref50/uniref50.fasta
cd data/datasets/blastdb/uniref50 && makeblastdb -in uniref50.fasta -parse_seqids -dbtype prot

echo "Extracting Uniref90..."
# gunzip -c data/raw/uniref/uniref90/uniref90.fasta.gz > data/datasets/blastdb/uniref90/uniref90.fasta
# cd data/datasets/blastdb/uniref90 && makeblastdb -in uniref90.fasta -parse_seqids -dbtype prot
