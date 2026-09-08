"""
Structural text renderings -> text embedding. 
env: hf
"""

from __future__ import annotations

from chemembed.registry import TEXT_MODELS
from chemembed.text.descriptions import load_descriptions
from generators._base import finish, load_inputs, parse_spec
from generators._text_embed import embed_texts


def main():
    spec = parse_spec()                      # variant: one of STRUCT_FIELDS
    model = TEXT_MODELS[spec.model]
    df = load_inputs()

    descriptions = load_descriptions()       # {compound_id: {field: text}}
    ids = df["compound_id"].tolist()

    missing = [c for c in ids if c not in descriptions]
    if missing:
        raise KeyError(f"{len(missing)} compounds absent from descriptions.json "
                       f"(e.g. {missing[:3]}) - rerun build_descriptions()")

    texts = [descriptions[c].get(spec.variant, "") for c in ids]
    n_empty = sum(1 for t in texts if not t.strip())
    print(f"[{spec}] {len(texts)} texts, {n_empty} empty, "
          f"median {sorted(len(t) for t in texts)[len(texts) // 2]} chars")

    matrix = embed_texts(model, texts)

    finish(spec, ids, matrix, meta={
        "hf_id": model.hf_id,
        "pooling": model.pooling,
        "field": spec.variant,
        "n_empty": n_empty,
    })


if __name__ == "__main__":
    main()