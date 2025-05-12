from subpred.util import load_data


def get_sequence_dataset(
    datasets_path: str,
    organism_ids: set = None,
    swissprot_only: bool = False,
    max_sequence_evidence_code: int = 2,
    additional_proteins: set = None,
    remove_proteins_without_gene_names:bool = True
):
    df_uniprot = load_data("uniprot", folder_path=datasets_path)
    if swissprot_only:
        df_uniprot = df_uniprot[df_uniprot.reviewed]
    df_uniprot = df_uniprot[df_uniprot.protein_existence <= max_sequence_evidence_code]
    if organism_ids:
        # ability to add proteins from other organism, e.g. when studying a plant protein by expressing it in S. cerevisiae
        df_uniprot = df_uniprot[
            df_uniprot.organism_id.isin(organism_ids)
            | df_uniprot.index.isin(
                additional_proteins if additional_proteins else set()
            )
        ]
    if remove_proteins_without_gene_names:
        df_uniprot = df_uniprot[~df_uniprot.gene_names.isnull()]
    df_uniprot = df_uniprot[
        ["sequence", "reviewed", "protein_existence", "organism_id", "protein_names"]
    ].drop_duplicates()
    return df_uniprot
