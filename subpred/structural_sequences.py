"""
@author: adenger
"""

from subpred.util import load_data


def get_3Di_sequences(
    accessions,
    remove_unequal_len: bool = True,
    dataset_name: str = "3Di_alphafold4",
    folder_path: str = "../data/datasets",
):
    # accessions: list
    # remove_unequal_len: only return if they match the length of the aa sequence
    sequences_3Di = load_data(dataset_name=dataset_name, folder_path=folder_path)
    sequences_3Di = sequences_3Di[sequences_3Di.index.isin(accessions)]
    accessions_3Di = [accession for accession in accessions if accession in sequences_3Di.index]
    sequences_3Di = sequences_3Di.loc[accessions_3Di]  # try to preserve order
    if remove_unequal_len:
        sequences_3Di = sequences_3Di[sequences_3Di.len_matches_aa_sequence]
    return sequences_3Di.sequence3Di
