"""chemembed artifacts -> the {drug_name: vector} dict CellFlow expects.

CellFlow looks up conditions by the value in obs["perturbation"], i.e. the
SciPlex drug NAME. chemembed keys embeddings by InChIKey. 
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from backbone_eval import config as cfg
from chemembed.artifacts import artifact_dir, load_embedding
from chemembed.config import COMPOUNDS


def drug_reps(spec: str, drug_names: list[str]) -> dict[str, np.ndarray]:
    """
    {drug_name: embedding} plus a zero vector for the control token.
    cellflow uses a control token for their vehicle cells; im following that. 0 vector means no perturbation applied. 
    """
    table = pd.read_csv(COMPOUNDS / f"{cfg.DATASET}.csv")
    name_to_id = dict(zip(table["dataset_drug_name"].str.strip(), table["compound_id"]))

    names = [n for n in drug_names if n in name_to_id]
    missing = sorted(set(drug_names) - set(names))
    if missing:
        raise KeyError(f"{len(missing)} compounds absent from the compound table: "
                       f"{missing[:5]}")

    ids = [name_to_id[n] for n in names]
    unique_ids = list(dict.fromkeys(ids))
    E = load_embedding(spec, unique_ids)           
    position = {c: i for i, c in enumerate(unique_ids)}

    reps = {n: E[position[name_to_id[n]]].astype(np.float32) for n in names}
    reps[cfg.CONTROL_TOKEN] = np.zeros(E.shape[1], dtype=np.float32)
    return reps


def available(specs: list[str]) -> list[str]:
    return [s for s in specs if (artifact_dir(s) / "embeddings.npy").exists()]