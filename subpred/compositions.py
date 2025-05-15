"""
@author: adenger
"""

from collections import Counter
from itertools import product
import numpy as np
import pandas as pd
from joblib.parallel import delayed, Parallel, cpu_count

AMINO_ACIDS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]

# 3di: ACDEFGHIKLMNPQRSTVWYX


def __kmer_composition(sequence: str, k: int, kmers_dict: dict) -> pd.Series:
    counter = Counter(kmers_dict)
    counter.update([sequence[i : i + k] for i in range(len(sequence) - k + 1)])
    return pd.Series([x / (len(sequence) - k + 1) for x in counter.values()])


def __kmer_composition_batch(
    sequences: pd.Series, k: int, kmers_dict: dict
) -> pd.DataFrame:
    return sequences.apply(__kmer_composition, k=k, kmers_dict=kmers_dict)


def calculate_comp(
    sequences: pd.Series, k: int = 2, alphabet: list = AMINO_ACIDS, n_threads: int = 1
) -> pd.DataFrame:
    assert k > 0
    assert n_threads != 0 and n_threads
    n_chunks = n_threads if n_threads > 0 else cpu_count() + n_threads + 1
    n_chunks = min(sequences.shape[0], n_chunks)
    # sequences_chunks = np.array_split(sequences, indices_or_sections=n_chunks)
    sequences_chunk_idxs = np.array_split(
        np.arange(sequences.shape[0]), indices_or_sections=n_chunks
    )
    sequences_chunks = [
        sequences.iloc[sequences_chunk_idx]
        for sequences_chunk_idx in sequences_chunk_idxs
    ]

    kmers = ["".join(x) for x in product(alphabet, repeat=k)]
    kmers_dict = {kmer: 0 for kmer in kmers}

    chunk_results = Parallel(n_jobs=n_threads)(
        delayed(__kmer_composition_batch)(sequences_chunk, k, kmers_dict)
        for sequences_chunk in sequences_chunks
    )
    df_kmer_frequencies = pd.concat(chunk_results)
    feature_name = "AAC" if k == 1 else "PAAC" if k == 2 else f"KMER{k}"
    df_kmer_frequencies.columns = [f"{feature_name}__" + kmer for kmer in kmers]
    return df_kmer_frequencies


def calculate_aac(sequences: pd.Series, n_threads: int = -1) -> pd.DataFrame:
    return calculate_comp(sequences=sequences, k=1, n_threads=n_threads)


def calculate_paac(sequences: pd.Series, n_threads: int = -1) -> pd.DataFrame:
    return calculate_comp(sequences=sequences, k=2, n_threads=n_threads)
