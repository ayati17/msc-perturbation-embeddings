"""

SMILES / SELFIES sequence models that load through plain transformers. 
env: hf

  chemberta  DeepChem/ChemBERTa-77M-MTR   SMILES,  encoder, mean pooling
  selformer  HUBioDataLab/SELFormer       SELFIES, encoder, mean pooling
  chemgpt    ncfrey/ChemGPT-1.2B          SELFIES, decoder, last-token + left padding
"""

from __future__ import annotations
import os
from chemembed.registry import STRING_MODELS
from generators._base import finish, load_inputs, parse_spec
from generators._text_embed import embed_texts

UNK_TOLERANCE = 0.0          


def main():
    spec = parse_spec()
    model = STRING_MODELS[spec.model]
    df = load_inputs()

    seqs = _sequences(df["smiles_canonical"].tolist(), model.input)
    unk_rate = _unk_rate(model, seqs)

    lengths = sorted(len(s) for s in seqs)
    print(f"[{spec}] {len(seqs)} {model.input} sequences, "
          f"median {lengths[len(lengths) // 2]} chars, unk rate {unk_rate:.4f}")

    matrix = embed_texts(model, seqs)

    finish(spec, df["compound_id"].tolist(), matrix, meta={
        "hf_id": model.hf_id,
        "input": model.input,
        "pooling": model.pooling,
        "padding_side": model.padding_side,
        "unk_rate": round(unk_rate, 6),
    })


def _sequences(smiles: list[str], kind: str) -> list[str]:
    if kind == "smiles":
        return smiles
    if kind != "selfies":
        raise ValueError(f"unknown input kind {kind!r}")

    import selfies as sf

    out, failed = [], []
    for s in smiles:
        try:
            out.append(sf.encoder(s))
        except Exception as e:
            failed.append((s, type(e).__name__))
    if failed:
        raise RuntimeError(
            f"{len(failed)} SMILES couldn't be SELFIES-encoded, e.g. {failed[:3]}. "
        )
    return out


def _unk_rate(model, seqs: list[str]) -> float:
    from transformers import AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model.hf_id, token=os.environ.get("HF_TOKEN"))
    unk_id = tokenizer.unk_token_id
    if unk_id is None:
        return 0.0

    encoded = tokenizer(seqs, truncation=True, max_length=model.max_length)
    total = sum(len(ids) for ids in encoded["input_ids"])
    unks = sum(ids.count(unk_id) for ids in encoded["input_ids"])
    rate = unks / total if total else 0.0

    if rate > UNK_TOLERANCE:
        worst = max(zip(encoded["input_ids"], seqs),
                    key=lambda p: p[0].count(unk_id))[1]
        raise RuntimeError(
            f"{model.hf_id}: {unks}/{total} tokens are <unk> ({rate:.2%}). "
            f"Vocab mismatch - check the selfies version. Worst sequence: {worst[:120]}"
        )
    return rate


if __name__ == "__main__":
    main()