import requests
import re
from requests.adapters import HTTPAdapter, Retry
import argparse
import pandas as pd

# Original from https://www.uniprot.org/help/api_queries

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


# TODO create urls dynamically
# TODO support for stream download
# TODO support for gzipped download

# def __create_uniprot_url(
#     base_url="https://rest.uniprot.org/uniprotkb/stream",
#     compressed="true",
#     fields=[
#         "accession",
#         "id",
#         "gene_names",
#         "protein_name",
#         "organism_name",
#         "organism_id",
#         "keywordid",
#         "keyword",
#         "go_id",
#         "go",
#         "xref_tcdb",
#         "protein_existence",
#         "sequence",
#         "fragment",
#     ],
#     format="tsv",
#     query="* AND (reviewed:true)",
# ):
#     """Create URL for use with curl"""
#     params = {
#         "compressed": compressed,
#         "fields": ",".join(fields),
#         "format": format,
#         "query": query,
#     }

#     return f"{base_url}?{urlencode(params, quote_via=quote)}"


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = int(response.headers["x-total-results"])
        yield response, total
        batch_url = get_next_link(response.headers)


# test_url = "https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Cgene_names%2Cprotein_name%2Corganism_name%2Corganism_id%2Ckeyword%2Ckeywordid%2Cgo_id%2Cprotein_existence%2Cfragment%2Csequence%2Cprotein_families%2Cxref_tcdb&format=tsv&query=%28%2A%29%20AND%20%28proteins_with%3A4%29%20AND%20%28model_organism%3A9606%29&size=500"
# url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=Insulin%20AND%20%28reviewed%3Atrue%29&size=500'


def download_dataset(url: str, output_path: str):
    df = pd.DataFrame()
    count = 0
    for batch, total in get_batch(url):
        print(
            f"downloading entries {count}-{min(total, count+500)} of {total} ({round(count / total * 100, 2)}%)...",
            end="\r",
        )
        count += 500
        lines = batch.text.splitlines()
        header = lines[0].split("\t")
        rows = []
        for line in batch.text.splitlines()[1:]:
            values = line.split("\t")
            rows.append(values)
        df_batch = pd.DataFrame(data=rows, columns=header)
        df_batch.set_index("Entry", inplace=True)
        df = pd.concat([df, df_batch])

    print(f"writing to file {output_path}...")
    print("done.")
    df.to_csv(output_path, sep="\t")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="Uniprot downloader",
        description="Download custom Uniprot dataset using pagination",
    )

    parser.add_argument("url", type=str)
    parser.add_argument("output_path", type=str)

    args = parser.parse_args()

    download_dataset(args.url, args.output_path)
