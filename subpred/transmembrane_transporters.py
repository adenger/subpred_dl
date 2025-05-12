from subpred.protein_dataset import get_sequence_dataset
from subpred.go_annotations import get_go_annotations_subset
from subpred.cdhit import cd_hit
# from subpred.chebi_annotations import get_go_chebi_annotations
import numpy as np
import pandas as pd
from subpred.util import load_data
import multiprocessing

# "human": 		9606
# "athaliana":	3702
# "ecoli": 		83333
# "yeast": 		559292


def get_transmembrane_transporter_dataset(
    organism_ids: set = None,
    swissprot_only: bool = False,
    datasets_path: str = "../data/datasets/",
    exclude_iea_go_terms: bool = False,
    max_sequence_evidence_code: int = 1,
    additional_proteins: set = None,
    # anatomical_entities_whitelist: set = None,
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
        anatomical_entities_whitelist (set, optional):
            Only keep proteins annotated with these CC GO terms. Defaults to None.
        remove_proteins_without_gene_names (bool, optional):
            Remove proteins that have not been annotated with a gene name in Uniprot. Defaults to True.

    Returns:
        df_sequences: Proteins dataset
        df_uniprot_goa: GO annotations for proteins
        df_go_chebi: ChEBI annotations for GO terms
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
        root_go_term="transmembrane transporter activity",
        inner_go_relations={"is_a"},
        namespaces_keep={"molecular_function"},
        proteins_subset=set(df_sequences.index),
        go_protein_qualifiers_filter_set={"enables"},
        annotations_evidence_codes_remove={"IEA"} if exclude_iea_go_terms else None,
    )
    # Filter sequences for those with transporter go annotations
    df_sequences = df_sequences[df_sequences.index.isin(df_uniprot_goa.Uniprot)]

    # Filter for cellular components
    # if anatomical_entities_whitelist:
    #     df_uniprot_goa_anatomical_entity = get_go_annotations_subset(
    #         datasets_path=datasets_path,
    #         root_go_term="cellular anatomical entity",
    #         inner_go_relations={"is_a"},
    #         namespaces_keep={"cellular_component"},
    #         proteins_subset=set(df_sequences.index),
    #         go_protein_qualifiers_filter_set={"located_in", "is_active_in"},
    #         annotations_evidence_codes_remove={"IEA"} if exclude_iea_go_terms else None,
    #     )
    #     anatomical_entities_protein_subset = df_uniprot_goa_anatomical_entity[
    #         df_uniprot_goa_anatomical_entity.go_term_ancestor.isin(
    #             anatomical_entities_whitelist
    #         )
    #     ].Uniprot.unique()

    #     df_sequences = df_sequences[
    #         df_sequences.index.isin(anatomical_entities_protein_subset)
    #     ]

    # Get chebi terms associated with go terms. Get them for ancestors, since go_id is subset of that.
    # df_go_chebi = get_go_chebi_annotations(
    #     dataset_path=datasets_path,
    #     go_ids_subset=set(df_uniprot_goa.go_id_ancestor),
    #     go_chebi_relations_subset={"has_primary_input", "has_participant"},
    #     filter_by_3star=False,
    #     add_ancestors=True,
    # )
    return df_sequences, df_uniprot_goa#, df_go_chebi


def get_stats(df_sequences, df_uniprot_goa, dataset_path="../data/datasets/"):

    df_sequences_merge = df_sequences.join(load_data("uniprot", folder_path=dataset_path)["gene_names"], how="left")
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
