"""
@author: adenger
"""

from subpred.cdhit import cd_hit
import pandas as pd
from subpred.structural_sequences import get_3Di_sequences
from subpred.cdhit import cd_hit


def cluster_sequences(dataset_full, identity_threshold: int):
    # dataset_full: tuple of df_sequences, df_uniprot_goa created with protein_go_datasets.py
    df_sequences, df_uniprot_goa = dataset_full[0].copy(), dataset_full[1].copy()
    cluster_representatives = set(
        cd_hit(df_sequences.sequence, identity_threshold=identity_threshold)
    )
    df_sequences = df_sequences[df_sequences.index.isin(cluster_representatives)]
    df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.index.isin(cluster_representatives)]
    return df_sequences, df_uniprot_goa


def filter_no_3Di_available(dataset_full, remove_unequal_len:bool=True):
    # dataset_full: tuple of df_sequences, df_uniprot_goa created with protein_go_datasets.py
    # remove_unequal_len: remove proteins where lengths of 3Di and AA sequence do not match
    # function to remove proteins that don't have all features available
    accessions_3Di_available = get_3Di_sequences(dataset_full[0].index, remove_unequal_len=remove_unequal_len).index
    dataset_full_filtered = (
        dataset_full[0][dataset_full[0].index.isin(accessions_3Di_available)].copy(),
        dataset_full[1][dataset_full[1].index.isin(accessions_3Di_available)].copy(),
    )
    return dataset_full_filtered


def get_classification_subset(dataset_full, go_terms: list):
    # dataset_full: tuple of df_sequences, df_uniprot_goa created with protein_go_datasets.py
    df_sequences, df_uniprot_goa = dataset_full[0].copy(), dataset_full[1].copy()
    df_uniprot_goa = (
        df_uniprot_goa[df_uniprot_goa.go_term_ancestor.isin(go_terms)][
            ["Uniprot", "go_term_ancestor"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    ).set_index("Uniprot")
    # Transporters that are annotated with both terms (might happen for other datasets)
    assert (
        not df_uniprot_goa.index.duplicated().any()
    ), "Some proteins are annotated with two or more of the GO terms simultaneously"
    df_sequences = df_sequences[
        df_sequences.index.isin(df_uniprot_goa.index)
    ]
    return df_sequences, df_uniprot_goa


def get_proteome_classification_subset(dataset_full, go_term: str):
    # dataset_full: tuple of df_sequences, df_uniprot_goa created with protein_go_datasets.py
    # go_term = "sugar transmembrane transporter activity"
    # go_term = "carbohydrate transmembrane transporter activity"
    df_sequences, df_uniprot_goa = dataset_full[0].copy(), dataset_full[1].copy()
    go_term_proteins = df_uniprot_goa[
        df_uniprot_goa.go_term_ancestor == go_term
    ].Uniprot.unique()
    background_proteins = df_uniprot_goa[
        ~df_uniprot_goa.go_term_ancestor.isin(go_term_proteins)
    ].Uniprot.unique()
    go_term_proteins_labels = [go_term] * len(go_term_proteins)
    background_proteins_labels = ["NOT|" + go_term] * len(background_proteins)

    df_task = pd.concat(
        [
            pd.DataFrame(
                data=go_term_proteins_labels,
                index=go_term_proteins,
                columns=["go_term_ancestor"],
            ),
            pd.DataFrame(
                data=background_proteins_labels,
                index=background_proteins,
                columns=["go_term_ancestor"],
            ),
        ],
        axis=0,
    )

    df_task.index = df_task.index.rename("Uniprot")
    df_sequences = df_sequences.loc[df_task.index]
    return df_sequences, df_task



