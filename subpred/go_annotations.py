"""
@author: adenger
"""

from subpred.util import load_data
import networkx as nx
import pandas as pd

EVIDENCE_CODE_TO_DESCRIPTION = {
    "IMP": "experimental_evidence",
    "IPI": "experimental_evidence",
    "IEP": "experimental_evidence",
    "IDA": "experimental_evidence",
    "EXP": "experimental_evidence",
    "IGI": "experimental_evidence",
    "HDA": "experimental_evidence_high_throughput",
    "HMP": "experimental_evidence_high_throughput",
    "HTP": "experimental_evidence_high_throughput",
    "HGI": "experimental_evidence_high_throughput",
    "HEP": "experimental_evidence_high_throughput",
    "IBA": "phylogenetically_inferred",
    "IBD": "phylogenetically_inferred",
    "IKR": "phylogenetically_inferred",
    "IRD": "phylogenetically_inferred",
    "ISS": "computational_analysis",
    "ISO": "computational_analysis",
    "ISA": "computational_analysis",
    "ISM": "computational_analysis",
    "IGC": "computational_analysis",
    "RCA": "computational_analysis",
    "NAS": "author_statement",
    "TAS": "author_statement",
    "IC": "curator_statement",
    "ND": "curator_statement",
    "IEA": "electronic_annotation",
}


def get_go_subgraph(
    graph_go,
    root_node: str = "GO:0022857",
    keys: set = {"is_a"},
    namespaces: set = {"molecular_function"},
):
    graph_go_subgraph = graph_go.subgraph(
        {
            node
            for node, namespace in graph_go.nodes(data="namespace")
            if namespace in namespaces
        }
    )
    graph_go_subgraph = graph_go_subgraph.subgraph(
        nx.ancestors(graph_go, root_node) | {root_node}
    )

    graph_go_subgraph = graph_go_subgraph.edge_subgraph(
        {
            (go1, go2, key)
            for go1, go2, key in graph_go_subgraph.edges(keys=True)
            if key in keys
        }
    )

    return graph_go_subgraph


def __get_id_update_dict(graph_go):
    go_id_update_dict = dict()
    for go_term, alt_ids in graph_go.nodes(data="alt_id"):
        if not alt_ids:
            go_id_update_dict[go_term] = go_term
            continue
        for alt_id in alt_ids:
            go_id_update_dict[alt_id] = go_term
    for go_term in graph_go.nodes():
        go_id_update_dict[go_term] = go_term
    return go_id_update_dict


def __get_go_annotations(
    datasets_path,
    go_id_update_dict,
    proteins_subset: set = None,
    go_ids_subset: set = None,
    qualifiers_keep: set = None,
    aspects_keep: set = None,
    evidence_codes_remove: set = None,
):
    df_uniprot_goa = load_data("goa", folder_path=datasets_path)

    # update go identifiers in annotation dataset to match go graph
    df_uniprot_goa["go_id"] = df_uniprot_goa.go_id.map(go_id_update_dict)
    df_uniprot_goa = df_uniprot_goa[~df_uniprot_goa.go_id.isnull()].reset_index(
        drop=True
    )

    # filtering out "not" annotations explicitly
    df_uniprot_goa = df_uniprot_goa[~df_uniprot_goa.qualifier.str.startswith("NOT")]

    # filtering for parameters
    if qualifiers_keep:
        df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.qualifier.isin(qualifiers_keep)]
    if aspects_keep:
        df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.aspect.isin(aspects_keep)]
    if evidence_codes_remove:
        df_uniprot_goa = df_uniprot_goa[
            ~df_uniprot_goa.evidence_code.isin(evidence_codes_remove)
        ]

    if proteins_subset:
        df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.Uniprot.isin(proteins_subset)]

    # filter annotations by protein subset
    if go_ids_subset:
        df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.go_id.isin(go_ids_subset)]

    # cleanup
    df_uniprot_goa = df_uniprot_goa.drop_duplicates().reset_index(drop=True)

    return df_uniprot_goa


def __add_ancestors(df_uniprot_goa, graph_go):
    df_uniprot_goa["go_id_ancestor"] = df_uniprot_goa.go_id.map(
        lambda go_id: {go_id} | set(nx.descendants(graph_go, go_id))
    )  # TODO this is very slow for larger graphs. multithreading did not improve running time by much, maybe try nx-cugraph

    df_uniprot_goa = df_uniprot_goa.explode("go_id_ancestor")
    df_uniprot_goa = df_uniprot_goa.reset_index(drop=True)

    return df_uniprot_goa


def get_go_annotations_subset(
    datasets_path: str,
    root_go_term: str,
    inner_go_relations: set = {"is_a"},
    namespaces_keep: set = {
        "biological_process",
        "molecular_function",
        "cellular_component",
    },
    proteins_subset: set = None,
    go_protein_qualifiers_filter_set: set = None,
    annotations_evidence_codes_remove: set = None,
) -> pd.DataFrame:
    """Creates go subset with protein annotations and ancestors of annotated terms.
    If a protein is annotated with "sodium ion uniporter activity",
    and  "sodium ion uniporter activity" is_a "monoatomic ion transmembrane transporter activity"
    according to GO, then the protein is also annotated with the latter. This increases the number of samples
    when training models with the GO terms as labels.

    Args:
        datasets_path (str):
            location of the pickles generated by preprocessing (data/datasets)
        root_go_term (str):
            only go terms that are children or ancestors of this term are included
        inner_go_relations (set):
            valid relations (edge labels) for searching for ancestors of root term in the go graph
        namespaces_keep (set):
            limit go graph to namespaces.
            options are: "biological_process""molecular_function", "cellular_component"
        proteins_subset (set, optional):
            subset the protein annotations to set of Uniprot accessions
            e.g.  to limit annotations to an organism
            Defaults to None, which means no filtering.
        go_protein_qualifiers_filter_set (set, optional):
            filter the valid interactions between proteins and go terms.
            for example "enables" or "located_in" or "part_of"
            Defaults to None, which means no filtering.
        annotations_evidence_codes_remove (set, optional):
            Filter out evidence codes, most commonly IEA. Defaults to None, which means no filtering.

    Returns:
        pd.DataFrame: Subset of go annotations, with added ancestors
    """

    # Example:

    # df_uniprot_goa = get_go_annotations_subset(
    #     datasets_path="../data/datasets/",
    #     root_go_term="transmembrane transporter activity",
    #     inner_go_relations={"is_a"},
    #     namespaces_keep={"molecular_function"},
    #     proteins_subset=set(df_sequences.index),
    #     go_protein_qualifiers_filter_set={"enables"},
    #     annotations_evidence_codes_remove={"IEA"},
    # )

    graph_go = load_data("go_obo", folder_path=datasets_path)
    go_id_to_term = {go_id: go_term for go_id, go_term in graph_go.nodes(data="name")}
    go_term_to_id = {go_term: go_id for go_id, go_term in graph_go.nodes(data="name")}
    go_id_update_dict = __get_id_update_dict(graph_go=graph_go)
    # This graph is not filtered by the protein set!
    graph_go_subgraph = get_go_subgraph(
        graph_go=graph_go,
        root_node=go_term_to_id[root_go_term],
        keys=inner_go_relations,
        namespaces=namespaces_keep,
    )
    namespace_to_aspect = {
        "biological_process": "P",
        "molecular_function": "F",
        "cellular_component": "C",
    }

    # get all go annotations, filtered by parameters, with updated ids
    df_uniprot_goa = __get_go_annotations(
        datasets_path=datasets_path,
        go_id_update_dict=go_id_update_dict,
        proteins_subset=proteins_subset,
        go_ids_subset=set(graph_go_subgraph.nodes()),
        qualifiers_keep=go_protein_qualifiers_filter_set,
        aspects_keep={namespace_to_aspect[namespace] for namespace in namespaces_keep},
        evidence_codes_remove=annotations_evidence_codes_remove,
    )
    # add ancestors
    df_uniprot_goa = __add_ancestors(
        df_uniprot_goa=df_uniprot_goa, graph_go=graph_go_subgraph
    )
    # add go terms
    df_uniprot_goa["go_term"] = df_uniprot_goa.go_id.map(go_id_to_term)
    df_uniprot_goa["go_term_ancestor"] = df_uniprot_goa.go_id_ancestor.map(
        go_id_to_term
    )
    # sort columns
    df_uniprot_goa = df_uniprot_goa[
        [
            "Uniprot",
            "qualifier",
            "go_id",
            "go_term",
            "evidence_code",
            "aspect",
            "go_id_ancestor",
            "go_term_ancestor",
        ]
    ]

    return df_uniprot_goa
