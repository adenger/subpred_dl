"""
@author: adenger
"""
from subpred.fasta import read_fasta, write_fasta
import tempfile
import subprocess
import os
import pandas as pd
import re



def __parse_cluster_file(file_path: str) -> pd.DataFrame:
    accession_pattern = re.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )
    cluster_file_header_pattern = re.compile(">Cluster ([0-9]+)\n")
    cluster_file_row_pattern = re.compile(
        "[0-9]+\t([0-9]+)aa, >([0-9A-Za-z]+)\.\.\..* ([0-9]{2,3}\.[0-9]{2}|\*).*\n"
    )
    records = []
    current_cluster = -1
    with open(file_path) as clstr_file:
        for row_str in clstr_file:
            if row_str.startswith(">"):
                current_cluster = int(
                    re.search(cluster_file_header_pattern, row_str).group(1)
                )
            else:
                matches = re.search(cluster_file_row_pattern, row_str)
                if not matches:
                    print(row_str)
                n_amino_acids, accession, percentage = (
                    matches.group(i) for i in [1, 2, 3]
                )
                if percentage == "*":
                    percentage = 100.00
                percentage = float(percentage)
                n_amino_acids = int(n_amino_acids)
                assert re.fullmatch(accession_pattern, accession), accession
                assert percentage != 0.0, row_str
                records.append([current_cluster, accession, percentage, n_amino_acids])
    df_clusters = pd.DataFrame.from_records(
        records,
        columns=["cluster", "accession", "identity_to_representative", "n_amino_acids"],
    )
    return df_clusters


def __auto_word_length(identity_threshold: int):
    if identity_threshold >= 70:
        return 5
    if identity_threshold >= 60:
        return 4
    if identity_threshold >= 50:
        return 3
    if identity_threshold >= 40:
        return 2
    raise ValueError("Invalid identity threshold: ", identity_threshold)


def __flatten_kwargs(**kwargs):
    kwargs_list = list()
    for k, v in kwargs.items():
        kwargs_list.append(k)
        kwargs_list.append(v)
    return kwargs_list


def cd_hit(
    sequences: pd.Series,
    identity_threshold: int,
    verbose: bool = True,
    executable_location: str = "cd-hit",
    n_threads: int = 4,
    memory: int = 4096,
    return_cluster_file: bool = False,
    **kwargs,
):
    fasta_data = list(
        zip([">" + ac for ac in sequences.index.tolist()], sequences.values.tolist())
    )

    with (
        tempfile.NamedTemporaryFile(suffix=".fasta") as tmp_fasta_in,
        tempfile.NamedTemporaryFile(suffix=".fasta") as tmp_fasta_out,
    ):
        write_fasta(fasta_file_name=tmp_fasta_in.name, fasta_data=fasta_data)

        word_length = __auto_word_length(identity_threshold)
        execution = [
            executable_location,
            "-i",
            tmp_fasta_in.name,
            "-o",
            tmp_fasta_out.name,
            "-c",
            str(identity_threshold / 100.0),
            "-n",
            str(word_length),
            "-T",
            str(n_threads),
            "-M",
            str(memory),
        ] + __flatten_kwargs(**kwargs)
        result = subprocess.run(
            execution, check=True, stdout=subprocess.PIPE, universal_newlines=True
        )
        if verbose:
            for line in result.stdout.split("\n"):
                if "finished" in line:
                    line = line.split()
                    print(
                        f"cd-hit: clustered {line[0]} sequences into {line[2]} clusters at threshold {identity_threshold}"
                    )
                    break

        if return_cluster_file:
            cluster_file_path = tmp_fasta_out.name + ".clstr"
            df_clusters = __parse_cluster_file(cluster_file_path)
            os.remove(tmp_fasta_out.name + ".clstr")
            return df_clusters
        else:
            fasta_data_clustered = read_fasta(tmp_fasta_out.name)
            os.remove(tmp_fasta_out.name + ".clstr")
            return [ac[1:] for ac, _ in fasta_data_clustered]
