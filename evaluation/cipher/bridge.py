#evaluation/cipher/bridge.py — feed chemembed artifacts to CIPHER pipeline 

from __future__ import annotations

import numpy as np
import pandas as pd

from chemembed.artifacts import artifact_dir, load_embedding
from chemembed.config import COMPOUNDS


def load_from_artifacts(specs: list[str], perturbation_names: np.ndarray,
                        dataset: str = "sciplex3") -> dict[str, np.ndarray]:
    table = pd.read_csv(COMPOUNDS / f"{dataset}.csv")
    name_to_id = dict(zip(table["dataset_drug_name"].str.strip(),
                          table["compound_id"]))

    names = pd.Series(perturbation_names).astype(str).str.strip()
    missing = sorted(set(names[~names.isin(name_to_id)]))
    if missing:
        raise KeyError(f"{len(missing)} compounds absent from {dataset}.csv: "
                       f"{missing[:5]}")

    ids = [name_to_id[n] for n in names]
    unique_ids = list(dict.fromkeys(ids))
    position = {c: i for i, c in enumerate(unique_ids)}
    take = np.array([position[c] for c in ids])

    out = {}
    for spec in specs:
        if not (artifact_dir(spec) / "embeddings.npy").exists():
            raise FileNotFoundError(f"no artifact for {spec}")
        E = load_embedding(spec, unique_ids)          # raises on anything absent
        out[spec.replace("/", "__")] = E[take].astype(np.float32)
    return out