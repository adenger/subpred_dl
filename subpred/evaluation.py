"""
@author: adenger
"""

def get_classification_subset(dataset_full, go_terms: list):
    df_sequences, df_uniprot_goa = dataset_full[0].copy(), dataset_full[1].copy()
    df_uniprot_goa = (
        df_uniprot_goa[df_uniprot_goa.go_term_ancestor.isin(go_terms)][
            ["Uniprot", "go_term_ancestor"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    ).set_index("Uniprot")
    # Transporters that are annotated with both terms (might happen for other datasets)
    assert not df_uniprot_goa.index.duplicated().any()
    df_sequences = df_sequences[df_sequences.index.isin(df_uniprot_goa.index)]
    return df_sequences, df_uniprot_goa

from subpred.cdhit import cd_hit

def cluster_sequences(dataset_full, identity_threshold: int):
    df_sequences, df_uniprot_goa = dataset_full[0].copy(), dataset_full[1].copy()
    cluster_representatives = set(
        cd_hit(df_sequences.sequence, identity_threshold=identity_threshold)
    )
    df_sequences = df_sequences[df_sequences.index.isin(cluster_representatives)]
    df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.index.isin(cluster_representatives)]
    return df_sequences, df_uniprot_goa