"""
@author: adenger
"""

from subpred.compositions import calculate_aac, calculate_paac
from subpred.pssm import calculate_pssm_feature
import pandas as pd
from sklearn.preprocessing import scale


def calculate_features(
    series_sequences: pd.Series,
    standardize_samples: bool = False,
    pssm_folder: str = "../data/datasets/pssm/",
    blastdb_folder: str = "../data/datasets/blastdb/",
    verbosity_pssm:int=1
):
    df_aac = calculate_aac(series_sequences)
    df_paac = calculate_paac(series_sequences)
    df_pssm_50_1 = calculate_pssm_feature(
        series_sequences,
        tmp_folder=pssm_folder + "pssm_uniref50_1it",
        blast_db=blastdb_folder + "uniref50/uniref50.fasta",
        iterations=1,
        psiblast_threads=-1,
        verbosity=verbosity_pssm,
        feature_name="PSSM_50_1",
    )
    df_pssm_50_3 = calculate_pssm_feature(
        series_sequences,
        tmp_folder=pssm_folder + "pssm_uniref50_3it",
        blast_db=blastdb_folder + "uniref50/uniref50.fasta",
        iterations=3,
        psiblast_threads=-1,
        verbosity=verbosity_pssm,
        feature_name="PSSM_50_3",
    )
    df_pssm_90_1 = calculate_pssm_feature(
        series_sequences,
        tmp_folder=pssm_folder + "pssm_uniref90_3it",
        blast_db=blastdb_folder + "uniref90/uniref90.fasta",
        iterations=1,
        psiblast_threads=-1,
        verbosity=verbosity_pssm,
        feature_name="PSSM_90_1",
    )
    df_pssm_90_3 = calculate_pssm_feature(
        series_sequences,
        tmp_folder=pssm_folder + "pssm_uniref90_3it",
        blast_db=blastdb_folder + "uniref90/uniref90.fasta",
        iterations=3,
        psiblast_threads=-1,
        verbosity=verbosity_pssm,
        feature_name="PSSM_90_3",
    )
    features_list = [
        df_aac,
        df_paac,
        df_pssm_50_1,
        df_pssm_50_3,
        df_pssm_90_1,
        df_pssm_90_3,
    ]
    if standardize_samples:
        features_list = [
            pd.DataFrame(
                data=scale(feature, axis=1),
                index=feature.index,
                columns=feature.columns,
            )
            for feature in features_list
        ]
    df_features = pd.concat(
        features_list,
        axis=1,
    )
    return df_features
