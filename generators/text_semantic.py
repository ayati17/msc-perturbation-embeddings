"""
Semantic text -> embedding. 
env: hf

Text comes from chemembed.text.semantic
notes for later: 
PubChem's description endpoint gave 0% DrugBank and 79% ChEBI, & the ChEBI text was 
structural restatement rather than pharmacology soo fall back on chemBL instead 

empty text gets embedded as the hf model's null vector 
count recorded in meta.json 
"""

from __future__ import annotations

from chemembed.registry import TEXT_MODELS
from chemembed.text.semantic import load_semantic, load_sources
from generators._base import finish, load_inputs, parse_spec
from generators._text_embed import embed_texts


def main():
    spec = parse_spec()                      # variant: "chembl" | "curated"
    model = TEXT_MODELS[spec.model]
    df = load_inputs()

    semantic = load_semantic(spec.variant)   # {compound_id: text}
    sources = load_sources()                 # {compound_id: "chembl"|"annotation"|"none"}

    ids = df["compound_id"].tolist()
    texts = [semantic.get(c, "") for c in ids]
    used = [sources.get(c, "none") for c in ids]

    n_empty = sum(1 for t in texts if not t.strip(". "))
    print(f"[{spec}] {len(texts)} texts, {n_empty} empty")

    matrix = embed_texts(model, texts)

    finish(spec, ids, matrix, meta={
        "hf_id": model.hf_id,
        "pooling": model.pooling,
        "variant": spec.variant,
        "n_empty": n_empty,
        "source_counts": {s: used.count(s) for s in sorted(set(used))},
    })


if __name__ == "__main__":
    main()