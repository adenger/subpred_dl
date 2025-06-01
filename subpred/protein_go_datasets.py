"""
@author: adenger
"""

from subpred.protein_dataset import get_sequence_dataset
from subpred.go_annotations import get_go_annotations_subset
from subpred.cdhit import cd_hit
import numpy as np
import pandas as pd
from subpred.util import load_data
import multiprocessing

# "human": 		9606
# "athaliana":	3702
# "ecoli": 		83333
# "yeast": 		559292

# TODO turn methods into dataclass
def get_uniprot_go_dataset(
    organism_ids: set = None,
    swissprot_only: bool = False,
    datasets_path: str = "../data/datasets/",
    max_sequence_evidence_code: int = 1,
    additional_proteins: set = None,
    remove_proteins_without_gene_names: bool = True,
    root_go_term: str = "transmembrane transporter activity",
    inner_go_relations: set = {"is_a"},
    go_uniprot_relations: set = {"enables"},
    namespaces_keep: set = {"molecular_function"},
    annotations_evidence_codes_remove: set = {"IEA"},
):
    """Creates the data that is used by many of the remaining methods as input

    Args:
        organism_ids (set, optional): Set of taxonomy identifiers, e.g. 9606 for human.
            No filtering occurs when set to None. Defaults to None.
        swissprot_only (bool, optional):
            Only keep proteins that have been manually curated by SwissProt. Defaults to False.
        datasets_path (str, optional):
            Path to the pickles from the preprocessing notebook. Defaults to "../data/datasets/".
        exclude_iea_go_terms (bool, optional):
            Whether to exclude GO term annotations that were inferred electronically. Defaults to False.
        max_sequence_evidence_code (int, optional):
            1 for evidence at protein level
            2 for evidence at protein or transcript level.
            Defaults to 1.
        additional_proteins (set, optional):
            Proteins to add to the dataset, for example from other organisms. Defaults to None.
        remove_proteins_without_gene_names (bool, optional):
            Remove proteins that have not been annotated with a gene name in Uniprot. Defaults to True.
        root_go_term (str, optional):
            GO term that is used as the root of the new GO subset/slim, e.g. "transmembrane transporter activity" or "plasma membrane".
            More abstract terms increase the size of the GO graph exponentially, and therefore the running time when searching for ancestors.
            If large subsets are required then optimize the get_go_annotations_subset.__add_ancestors method, e.g. with CUDA or joblib
        inner_go_relations (set, optional):
            Accepted relations between GO terms, such as "is_a" or "part_of", following the OBO format (https://geneontology.org/docs/ontology-relations/)
        go_uniprot_relations (set, optional):
            Accepted relations between GO terms and Uniprot proteins, such as "enables" or "located_in". Defined by the GAF 2.2 format.
        namespaces_keep (set, optional):
            Which aspect to keep GO terms from. Can contain "molecular_function", "cellular_component", and "biological_process"
        annotations_evidence_codes_remove (set, optional):
            Filter for GO-Uniprot annotations, according to the quality of evidence available. See https://geneontology.org/docs/guide-go-evidence-codes/


    Returns:
        df_sequences: Proteins dataset
        df_uniprot_goa: GO annotations for proteins
    """

    # First, get all sequences with filtering criteria:
    df_sequences = get_sequence_dataset(
        datasets_path=datasets_path,
        organism_ids=organism_ids,
        swissprot_only=swissprot_only,
        max_sequence_evidence_code=max_sequence_evidence_code,
        additional_proteins=additional_proteins,
        remove_proteins_without_gene_names=remove_proteins_without_gene_names,
    )

    # Get GO annotations from subset of transmembrane transporter go terms
    df_uniprot_goa = get_go_annotations_subset(
        datasets_path=datasets_path,
        root_go_term=root_go_term,
        inner_go_relations=inner_go_relations,
        namespaces_keep=namespaces_keep,
        proteins_subset=set(df_sequences.index),
        go_protein_qualifiers_filter_set=go_uniprot_relations,
        annotations_evidence_codes_remove=annotations_evidence_codes_remove,
    )
    # Filter sequences for those with transporter go annotations
    df_sequences = df_sequences[df_sequences.index.isin(df_uniprot_goa.Uniprot)]

    return df_sequences, df_uniprot_goa


def get_transmembrane_transporter_dataset(
    organism_ids: set = None,
    swissprot_only: bool = False,
    datasets_path: str = "../data/datasets/",
    exclude_iea_go_terms: bool = False,
    max_sequence_evidence_code: int = 1,
    additional_proteins: set = None,
    remove_proteins_without_gene_names: bool = True,
):
    """Creates the data that is used by many of the remaining methods as input

    Args:
        organism_ids (set, optional): Set of taxonomy identifiers, e.g. 9606 for human.
            No filtering occurs when set to None. Defaults to None.
        swissprot_only (bool, optional):
            Only keep proteins that have been manually curated by SwissProt. Defaults to False.
        datasets_path (str, optional):
            Path to the pickles from the preprocessing notebook. Defaults to "../data/datasets/".
        exclude_iea_go_terms (bool, optional):
            Whether to exclude GO term annotations that were inferred electronically. Defaults to False.
        max_sequence_evidence_code (int, optional):
            1 for evidence at protein level
            2 for evidence at protein or transcript level.
            Defaults to 1.
        additional_proteins (set, optional):
            Proteins to add to the dataset, for example from other organisms. Defaults to None.
        remove_proteins_without_gene_names (bool, optional):
            Remove proteins that have not been annotated with a gene name in Uniprot. Defaults to True.

    Returns:
        df_sequences: Proteins dataset
        df_uniprot_goa: GO annotations for proteins
    """
    df_sequences, df_uniprot_goa = get_uniprot_go_dataset(
        organism_ids=organism_ids,
        swissprot_only=swissprot_only,
        datasets_path=datasets_path,
        max_sequence_evidence_code=max_sequence_evidence_code,
        additional_proteins=additional_proteins,
        remove_proteins_without_gene_names=remove_proteins_without_gene_names,
        root_go_term="transmembrane transporter activity",
        inner_go_relations={"is_a"},
        go_uniprot_relations={"enables"},
        namespaces_keep={"molecular_function"},
        annotations_evidence_codes_remove={"IEA"} if exclude_iea_go_terms else None,
    )

    return df_sequences, df_uniprot_goa


def count_children(df_uniprot_goa, go_term: str, **kwargs):
    # df_goa = dataset_organism[1]
    matching_go_ids = df_uniprot_goa[df_uniprot_goa.go_term == go_term].go_id.drop_duplicates()
    assert len(matching_go_ids == 1), matching_go_ids
    go_id = str(matching_go_ids.iloc[0])

    network = load_data("go_obo", **kwargs)

    predecessors = set(network.predecessors(go_id)) | {go_id}

    return (
        df_uniprot_goa[df_uniprot_goa.go_id_ancestor.isin(predecessors)][
            ["Uniprot", "go_term_ancestor"]
        ]
        .drop_duplicates()
        .groupby("go_term_ancestor")
        .count()
        .sort_values("Uniprot", ascending=False)
    )

def get_interpro_annotations(df_uniprot_goa, type: str, **kwargs):
    assert type in [
        "Family",
        "Domain",
        "Repeat",
        "Homologous_superfamily",
        "Conserved_site",
        "Binding_site",
        "Active_site",
        "PTM",
    ]
    df_interpro = load_data("interpro", **kwargs)
    df_interpro = (
        df_interpro[df_interpro.Uniprot.isin(df_uniprot_goa.index.unique()) & (df_interpro.type == type)]
        .sort_values("Uniprot")
        .reset_index(drop=True)
    )
    return df_interpro


def get_stats(df_sequences, df_uniprot_goa, dataset_path="../data/datasets/"):

    df_sequences_merge = df_sequences.join(
        load_data("uniprot", folder_path=dataset_path)["gene_names"], how="left"
    )
    df_sequences_merge["has_gene_name"] = ~df_sequences_merge.gene_names.isnull()
    df_sequences_merge = df_sequences_merge.drop(
        ["gene_names", "protein_names"], axis=1
    )

    df_sequences_merge = df_sequences_merge.reset_index().drop_duplicates()
    df_sequences_merge

    df_sequences_goa_merged = pd.merge(
        df_sequences_merge[
            ["Uniprot", "reviewed", "protein_existence", "has_gene_name"]
        ],
        df_uniprot_goa[["Uniprot", "evidence_code", "go_term_ancestor"]],
        on="Uniprot",
        how="inner",
    )
    df_sequences_goa_merged
    df_sequences_goa_merged["evidence_code"] = df_sequences_goa_merged[
        "evidence_code"
    ].transform(lambda x: "computational" if x == "IEA" else "experiment")
    df_sequences_goa_merged["protein_existence_evidence"] = df_sequences_goa_merged[
        "protein_existence"
    ].map({1: "protein_level", 2: "transcript_level"})
    df_sequences_goa_merged = df_sequences_goa_merged.drop("protein_existence", axis=1)
    df_sequences_goa_merged["clustering"] = "None"
    cdhit_cores = min(multiprocessing.cpu_count(), 12)

    for thresh in [50, 70, 90, 100]:
        cluster_representatives = cd_hit(
            df_sequences.sequence, identity_threshold=thresh, n_threads=cdhit_cores
        )

        df_sequences_goa_merged_clustered = (
            df_sequences_goa_merged[
                df_sequences_goa_merged.Uniprot.isin(cluster_representatives)
            ]
            .drop("clustering", axis=1)
            .assign(clustering=thresh)
            .drop_duplicates()
        )

        df_sequences_goa_merged = pd.concat(
            [df_sequences_goa_merged, df_sequences_goa_merged_clustered]
        ).reset_index(drop=True)
    df_sequences_goa_merged = df_sequences_goa_merged.rename(
        columns={
            "evidence_code": "go_evidence",
            "reviewed": "swissprot_reviewed",
            "go_term_ancestor": "go_term",
        }
    )

    df_sequences_goa_merged = df_sequences_goa_merged.drop_duplicates()

    df_stats_transporters = (
        df_sequences_goa_merged.drop("go_term", axis=1)
        .drop_duplicates()
        .groupby(
            [
                "swissprot_reviewed",
                "has_gene_name",
                "go_evidence",
                "protein_existence_evidence",
                "clustering",
            ]
        )
        .apply(np.unique)
        .apply(len)
        .to_frame("n_transporters")
    )

    df_stats_go = (
        df_sequences_goa_merged.drop("Uniprot", axis=1)
        .drop_duplicates()
        .groupby(
            [
                "swissprot_reviewed",
                "has_gene_name",
                "go_evidence",
                "protein_existence_evidence",
                "clustering",
            ]
        )
        .apply(np.unique)
        .apply(len)
        .to_frame("n_terms")
    )
    return df_stats_transporters.join(df_stats_go)
