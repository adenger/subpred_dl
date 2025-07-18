"""
@author: adenger
"""

# import os
from pathlib import Path
import pandas as pd
import networkx as nx
import pickle

# central methods for loading and saving pickles, since different data structures need different methods for that.


def save_data(
    dataset: pd.DataFrame | nx.Graph,
    dataset_name: str,
    folder_path: str = "../data/datasets",
    method: str = "auto",
    **kwargs,
):
    if method == "auto":
        if isinstance(dataset, pd.DataFrame):
            method = "pickle"
        elif isinstance(dataset, nx.Graph):
            method = "gpickle"
    match (method):
        case "pickle":
            dataset.to_pickle(f"{folder_path}/{dataset_name}.pickle", **kwargs)
        case "gpickle":
            with open(f"{folder_path}/{dataset_name}.gpickle", "wb") as pickle_file:
                pickle.dump(dataset, pickle_file)
        # case "pickle_gz":
        #     df.to_pickle(
        #         f"{folder_path}/{dataset_name}.pickle_gz",
        #         compression={"method": "gzip", "compresslevel": 1, "mtime": 1},
        #     )


def load_data(dataset_name: str, folder_path: str = "../data/datasets", **kwargs):
    for file_path in Path(folder_path).iterdir():
        file_name = file_path.stem
        if file_name != dataset_name:
            continue
        # file name is dataset name
        file_extension = file_path.suffix[1:]
        match (file_extension):
            case "pickle" | "pickle_gz":
                return pd.read_pickle(file_path, **kwargs)
            case "gpickle":
                with open(file_path, "rb") as pickle_file:
                    return pickle.load(pickle_file)


def save_results(
    dataset: pd.DataFrame | nx.Graph,
    dataset_name: str,
    folder_path: str = "../data/results",
    method: str = "auto",
    **kwargs,
):
    save_data(
        dataset=dataset,
        dataset_name=dataset_name,
        folder_path=folder_path,
        method=method,
        **kwargs,
    )


# import numpy as np
# from sklearn.pipeline import make_pipeline
# from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA
# from sklearn.feature_selection import f_classif, chi2
# from sklearn.cluster import AgglomerativeClustering
# Functions used my multiple notebooks and other functions

# def get_protein_aac_stats(df_aac, accession):
#     # compare the feature values of accession to the other features in the dataset
#     return pd.DataFrame(
#         data=[
#             [
#                 aa,
#                 df_aac.loc[accession, aa],
#                 df_aac.loc[df_aac.index != accession, aa].mean(),
#                 df_aac.loc[df_aac.index != accession, aa].std(),
#                 abs(
#                     df_aac.loc[df_aac.index != accession, aa].mean()
#                     - df_aac.loc[accession, aa]
#                 ),
#             ]
#             for aa in df_aac.columns
#         ],
#         columns=["AA", accession, "Mean", "Std", "Diff"],
#     ).sort_values("Diff", ascending=False).set_index("AA")


# def perform_pca(df_feature, labels: pd.Series, n_components: int = 2):
#     pipe = make_pipeline(StandardScaler(), PCA(n_components=n_components))
#     df_pca = pd.DataFrame(
#         data=pipe.fit_transform(df_feature),
#         columns=[f"PC{pc}" for pc in range(1, n_components + 1)],
#         index=df_feature.index,
#     )
#     df_pca[labels.name] = labels
#     return df_pca


# def get_feature_score(df_feature, labels: pd.Series, method: str = "f_classif", remove_zero_variance:bool = False):
#     if remove_zero_variance:
#         # does not overwrite original dataframe
#         df_feature = df_feature.loc[:,(df_feature.var() != 0)]
#     func = None
#     if method == "f_classif":
#         func = f_classif
#     elif method == "chi2":
#         func = chi2
#     else:
#         raise ValueError(f"Inalid method: {method}")

#     df_score = pd.DataFrame(
#         {
#             "Feature": df_feature.columns.tolist(),
#             "Normalized score": func(df_feature, labels)[0],
#             "Measure": [f"Feature importance ({method})"] * df_feature.shape[1],
#         }
#     )
#     df_score["Normalized score"] = df_score["Normalized score"] / sum(
#         df_score["Normalized score"]
#     )
#     return df_score

# def get_clusters(df_features, n_clusters=2):
#     # seems to be deterministic
#     cluster = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward")

#     cluster.fit(df_features)

#     cluster_list = []
#     for label in np.unique(cluster.labels_):
#         cluster_list.append(df_features.index[cluster.labels_ == label].tolist())

#     return cluster_list

# def get_protein_feature_stats(df_feature, accession):
#     accession_values = df_feature.loc[accession]
#     df_feature_not_accession = df_feature.loc[df_feature.index != accession]
#     df_stats_protein = pd.concat(
#         [
#             accession_values,
#             df_feature_not_accession.mean().rename("mean_feature"),
#             df_feature_not_accession.std().rename("std_feature"),
#         ],
#         axis=1,
#     )
#     df_stats_protein = df_stats_protein.assign(
#         diff=abs(df_stats_protein[accession] - df_stats_protein.mean_feature)
#     )
#     return df_stats_protein.sort_values("diff", ascending=False)
