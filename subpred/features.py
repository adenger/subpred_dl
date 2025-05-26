"""
@author: adenger
"""

from subpred.compositions import calculate_aac, calculate_paac
from subpred.pssm import calculate_pssm_feature
import pandas as pd
from sklearn.preprocessing import scale
from sklearn.preprocessing import LabelEncoder
from subpred.pssm import calculate_pssm_feature
from subpred.compositions import calculate_comp, ALPHABET_3DI, AMINO_ACIDS
from subpred.structural_sequences import get_3Di_sequences
from subpred.embeddings import get_nlp_features
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale


def get_ml_datasets(features_list: list, series_labels: pd.Series):
    # features_list: list of dataframes, where rows are samples and columns features
    ml_datasets = [
        (
            feature_name,
            *get_ml_dataset(df_features=df_features, series_labels=series_labels),
        )
        for feature_name, df_features in features_list
    ]
    return ml_datasets


def get_ml_dataset(df_features: pd.DataFrame, series_labels: pd.Series):
    # df_features, series_labels from features.get_features
    assert not series_labels.index.duplicated().any()
    assert (df_features.index == series_labels.index).all()

    sample_names = df_features.index.to_numpy()
    feature_names = df_features.columns.to_numpy()

    label_encoder = LabelEncoder()

    X = df_features.loc[sample_names].to_numpy()
    y_str = series_labels.loc[sample_names].to_numpy().ravel()
    y = label_encoder.fit_transform(y_str)

    return X, y, sample_names, feature_names


def get_features(dataset_full: tuple, include_pssm_features: bool = True):
    # dataset_full: generated with protein_go_datasets.py
    # Can take a long time if cache is empty
    df_sequences, df_uniprot_goa = dataset_full[0].copy(), dataset_full[1].copy()
    series_sequences = df_sequences.sequence
    series_accessions = df_sequences.index

    # 3Di sequences
    sequences_3Di = get_3Di_sequences(
        series_accessions,
        remove_unequal_len=True,
        dataset_name="3Di_alphafold4",
        folder_path="../data/datasets",
    )

    # Are there as many 3Di sequences as AA sequences? TODO If yes, maybe take intersection
    assert len(sequences_3Di) == len(series_sequences)
    assert (sequences_3Di.index == series_sequences.index).all()

    # original sequences features
    df_aac = calculate_comp(series_sequences, k=1, alphabet=AMINO_ACIDS)
    df_paac = calculate_comp(series_sequences, k=2, alphabet=AMINO_ACIDS)
    df_kmer3 = calculate_comp(series_sequences, k=3, alphabet=AMINO_ACIDS)

    df_3Di_AAC = calculate_comp(sequences=sequences_3Di, k=1, alphabet=ALPHABET_3DI)
    df_3Di_PAAC = calculate_comp(sequences=sequences_3Di, k=2, alphabet=ALPHABET_3DI)
    df_3Di_KMER3 = calculate_comp(sequences=sequences_3Di, k=3, alphabet=ALPHABET_3DI)

    # combining aa and 3di kmers
    print(df_3Di_AAC.index)
    print(df_aac.index)
    df_KMER1_COMBINED = pd.concat([df_3Di_AAC, df_aac], axis=1)
    df_KMER2_COMBINED = pd.concat([df_3Di_PAAC, df_paac], axis=1)
    df_KMER3_COMBINED = pd.concat([df_3Di_KMER3, df_kmer3], axis=1)

    # AA Embeddings
    df_embeddings_prott5_AA = get_nlp_features(
        sequences=series_sequences,
        model="protT5",
        sequence_type="AA",
        half_precision=True,
    )
    df_embeddings_prostt5_AA = get_nlp_features(
        sequences=series_sequences,
        model="prostT5",
        sequence_type="AA",
        half_precision=True,
    )
    # 3Di Embeddings
    df_embeddings_prott5_3Di = get_nlp_features(
        sequences=sequences_3Di,
        model="prostT5",
        sequence_type="3Di",
        half_precision=True,
    )

    np.random.seed(0)
    df_dummy_feature = pd.DataFrame(
        np.random.rand(len(series_accessions), 1024),
        index=series_accessions,
        columns=[f"dummy{i}" for i in range(1024)],
    )

    features_list = [
        ("DUMMY", df_dummy_feature),
        ("AAC", df_aac),
        ("PAAC", df_paac),
        ("AA_KMER3", df_kmer3),
        # ("PSSM_50_1", df_pssm_50_1),
        # ("PSSM_50_3", df_pssm_50_3),
        # ("PSSM_90_1", df_pssm_90_1),
        # ("PSSM_90_3", df_pssm_90_3),
        # ("PSSM_META", df_pssm_meta),
        # ("META", df_meta),
        # ("META_STD", df_meta_std),
        ("3Di_COMP", df_3Di_AAC),
        ("3Di_KMER2", df_3Di_PAAC),
        ("3Di_KMER3", df_3Di_KMER3),
        ("COMB_KMER1", df_KMER1_COMBINED),
        ("COMB_KMER2", df_KMER2_COMBINED),
        ("COMB_KMER3", df_KMER3_COMBINED),
        ("PROTT5_AA", df_embeddings_prott5_AA),
        ("PROSTT5_AA", df_embeddings_prostt5_AA),
        ("PROSTT5_3DI", df_embeddings_prott5_3Di),
    ]

    if include_pssm_features:
        pssm_folder = "../data/datasets/pssm/"
        blastdb_folder = "../data/datasets/blastdb/"
        verbosity_pssm = 1  # only print if no pssm found
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
        df_pssm_meta = pd.concat(
            [df_pssm_50_1, df_pssm_50_3, df_pssm_90_1, df_pssm_90_3], axis=1
        )

        # df_meta = calculate_features(series_sequences=series_sequences, standardize_samples=False)
        # df_meta_std = calculate_features(series_sequences=series_sequences, standardize_samples=True)
        # Meta feature from previous study
        df_meta = pd.concat(
            [
                df_aac,
                df_paac,
                df_pssm_50_1,
                df_pssm_50_3,
                df_pssm_90_1,
                df_pssm_90_3,
            ],
            axis=1,
        )
        # Meta fingerprint from previous study
        df_meta_std = pd.concat(
            [
                pd.DataFrame(
                    data=scale(feature.copy(), axis=1),
                    index=feature.index,
                    columns=feature.columns,
                )
                for feature in [
                    df_aac,
                    df_paac,
                    df_pssm_50_1,
                    df_pssm_50_3,
                    df_pssm_90_1,
                    df_pssm_90_3,
                ]
            ],
            axis=1,
        )

        features_list.extend(
            [
                # ("DUMMY", df_dummy_feature),
                # ("AAC", df_aac),
                # ("PAAC", df_paac),
                # ("AA_KMER3", df_kmer3),
                ("PSSM_50_1", df_pssm_50_1),
                ("PSSM_50_3", df_pssm_50_3),
                ("PSSM_90_1", df_pssm_90_1),
                ("PSSM_90_3", df_pssm_90_3),
                ("PSSM_META", df_pssm_meta),
                ("META", df_meta),
                ("META_STD", df_meta_std),
                # ("3Di_COMP", df_3Di_AAC),
                # ("3Di_KMER2", df_3Di_PAAC),
                # ("3Di_KMER3", df_3Di_KMER3),
                # ("COMB_KMER1", df_KMER1_COMBINED),
                # ("COMB_KMER2", df_KMER2_COMBINED),
                # ("COMB_KMER3", df_KMER3_COMBINED),
                # ("PROTT5_AA", df_embeddings_prott5_AA),
                # ("PROSTT5_AA", df_embeddings_prostt5_AA),
                # ("PROSTT5_3DI", df_embeddings_prott5_3Di),
            ]
        )
    features_list = [
        (feature_name, df_feature.loc[series_accessions])
        for feature_name, df_feature in features_list
    ]
    series_labels = df_uniprot_goa.loc[series_accessions].go_term_ancestor
    return features_list, series_labels


# def calculate_features(
#     series_sequences: pd.Series,
#     standardize_samples: bool = False,
#     pssm_folder: str = "../data/datasets/pssm/",
#     blastdb_folder: str = "../data/datasets/blastdb/",
#     verbosity_pssm:int=1
# ):
#     df_aac = calculate_aac(series_sequences)
#     df_paac = calculate_paac(series_sequences)
#     df_pssm_50_1 = calculate_pssm_feature(
#         series_sequences,
#         tmp_folder=pssm_folder + "pssm_uniref50_1it",
#         blast_db=blastdb_folder + "uniref50/uniref50.fasta",
#         iterations=1,
#         psiblast_threads=-1,
#         verbosity=verbosity_pssm,
#         feature_name="PSSM_50_1",
#     )
#     df_pssm_50_3 = calculate_pssm_feature(
#         series_sequences,
#         tmp_folder=pssm_folder + "pssm_uniref50_3it",
#         blast_db=blastdb_folder + "uniref50/uniref50.fasta",
#         iterations=3,
#         psiblast_threads=-1,
#         verbosity=verbosity_pssm,
#         feature_name="PSSM_50_3",
#     )
#     df_pssm_90_1 = calculate_pssm_feature(
#         series_sequences,
#         tmp_folder=pssm_folder + "pssm_uniref90_3it",
#         blast_db=blastdb_folder + "uniref90/uniref90.fasta",
#         iterations=1,
#         psiblast_threads=-1,
#         verbosity=verbosity_pssm,
#         feature_name="PSSM_90_1",
#     )
#     df_pssm_90_3 = calculate_pssm_feature(
#         series_sequences,
#         tmp_folder=pssm_folder + "pssm_uniref90_3it",
#         blast_db=blastdb_folder + "uniref90/uniref90.fasta",
#         iterations=3,
#         psiblast_threads=-1,
#         verbosity=verbosity_pssm,
#         feature_name="PSSM_90_3",
#     )
#     features_list = [
#         df_aac,
#         df_paac,
#         df_pssm_50_1,
#         df_pssm_50_3,
#         df_pssm_90_1,
#         df_pssm_90_3,
#     ]
#     if standardize_samples:
#         features_list = [
#             pd.DataFrame(
#                 data=scale(feature, axis=1),
#                 index=feature.index,
#                 columns=feature.columns,
#             )
#             for feature in features_list
#         ]
#     df_features = pd.concat(
#         features_list,
#         axis=1,
#     )
#     return df_features
