"""
@author: adenger
"""

from transformers import T5Tokenizer, T5EncoderModel
import torch
import pandas as pd
import numpy as np
from pathlib import Path


def generate_ProtT5_embeddings(sequences: pd.DataFrame, half_precision: bool = True):
    assert sequences.str.fullmatch("^[ACDEFGHIKLMNPQRSTVWY]+").all()
    weights_dtype = torch.float16 if half_precision else torch.float32

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    model_name = "Rostlab/prot_t5_xl_half_uniref50-enc"

    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_name, torch_dtype=weights_dtype).to(
        device
    )

    for accession, sequence in sequences.items():
        sequence_split = " ".join(sequence)
        ids = tokenizer(
            sequence_split,
            add_special_tokens=True,
            return_tensors="pt",
            padding="longest",
        ).to(device)
        with torch.no_grad():
            embedding_repr = model(ids.input_ids, attention_mask=ids.attention_mask)
        emb_single = embedding_repr.last_hidden_state[0, 0 : len(sequence)]
        emb_per_protein_single = emb_single.mean(dim=0).detach().cpu().numpy()
        torch.cuda.empty_cache()
        yield accession, emb_per_protein_single


def generate_ProstT5_embeddings(
    sequences: pd.DataFrame, sequence_type: str = "3Di", half_precision: bool = True
):
    weights_dtype = torch.float16 if half_precision else torch.float32
    model = "Rostlab/ProstT5_fp16"
    match sequence_type:
        case "3Di":
            special_token = "<fold2AA>"
            assert sequences.str.islower().all(), "3Di sequences need to be lower case"
            assert sequences.str.fullmatch("^[acdefghiklmnpqrstvwy]+").all()
        case "AA":
            special_token = "<AA2fold>"
            assert sequences.str.isupper().all(), "AA sequences need to be upper case"
            assert sequences.str.fullmatch("^[ACDEFGHIKLMNPQRSTVWY]+").all()
        case other:
            raise ValueError(f"unknown sequence type {other}")

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    tokenizer = T5Tokenizer.from_pretrained(model, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(
        model, torch_dtype=weights_dtype
    ).to(device)

    for accession, sequence in sequences.items():
        sequence_split = special_token + " " + " ".join(sequence)
        ids = tokenizer(
            sequence_split,
            add_special_tokens=True,
            return_tensors="pt",
            padding="longest",
        ).to(device)
        with torch.no_grad():
            embedding_repr = model(ids.input_ids, attention_mask=ids.attention_mask)
        emb_single = embedding_repr.last_hidden_state[0, 1 : len(sequence) + 1]
        emb_per_protein_single = emb_single.mean(dim=0).detach().cpu().numpy()
        torch.cuda.empty_cache()
        yield accession, emb_per_protein_single


def get_nlp_features(
    sequences: pd.Series,
    model: str = "protT5",
    sequence_type: str = "AA",
    half_precision: bool = True,
    cache_folder: Path = Path("../data/datasets/embeddings/"),
):
    # supports caching of results for faster loading
    # sequences: index uniprot accession, value sequence,
    #            lowercase for 3Di, uppercase for amino acid, only 20 AAs allowed
    # model: one of protT5, prostT5, TODO more?
    # sequence_type: AA or 3Di
    # half_precision: if false, weights and embeddings are fp32 instead of fp16,
    #                 uses more resources but more digits after comma. cpu only supports fp32
    assert cache_folder.exists() and cache_folder.is_dir()

    feature_name = f"{model}_{sequence_type}_{"fp16" if half_precision else "fp32"}"
    embedding_file_paths = sequences.index.to_series().map(
        lambda accession: cache_folder / f"{accession}_{feature_name}.npy"
    )

    sequences_noembeddings = sequences[~embedding_file_paths.map(lambda p: p.exists())]

    if len(sequences_noembeddings) > 0:
        match model:
            case "protT5":
                embeddings_res = generate_ProtT5_embeddings(
                    sequences=sequences_noembeddings, half_precision=half_precision
                )
            case "prostT5":
                embeddings_res = generate_ProstT5_embeddings(
                    sequences=sequences_noembeddings,
                    sequence_type=sequence_type,
                    half_precision=half_precision,
                )
                # TODO save
            case other:
                raise ValueError(f"Model not known: {other}")

        # save to cache
        for accession, embedding in embeddings_res:
            np.save(embedding_file_paths[accession], embedding)

        del embeddings_res

    torch.cuda.empty_cache()
    embeddings = np.array([np.load(x) for x in embedding_file_paths.values])
    df_embeddings = pd.DataFrame(data=embeddings, index=embedding_file_paths.index)
    df_embeddings.columns = [
        f"{feature_name}_{dim_count}" for dim_count in range(df_embeddings.shape[1])
    ]
    # prevent overflow warnings during printing, automatically select the best type
    df_embeddings = df_embeddings.convert_dtypes(convert_floating=True)

    return df_embeddings
