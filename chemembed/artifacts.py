"""
Artifact layout, one directory per spec:
data/artifacts/<spec>/embeddings.npy    (n_compounds, dim) float32, NOT NORMALIZED!!!!
                        index.parquet     compound_id, is_control, row order == npy
                        meta.json         
Vectors are stored w raw formats. Normalised at load time (load_embedding)
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from chemembed.config import ARTIFACTS
from chemembed.registry import Spec


def artifact_dir(spec: Spec | str) -> Path:
    return ARTIFACTS / str(spec)


def exists(spec: Spec | str) -> bool:
    d = artifact_dir(spec)
    return (d / "embeddings.npy").exists() and (d / "index.parquet").exists()


def write_artifact(
    spec: Spec | str,
    compound_ids: list[str],
    matrix: np.ndarray,
    *,
    is_control: list[bool] | None = None,
    meta: dict | None = None,
) -> Path:
    if len(compound_ids) != matrix.shape[0]:
        raise ValueError(f"{len(compound_ids)} ids vs {matrix.shape[0]} rows")
    if len(set(compound_ids)) != len(compound_ids):
        raise ValueError("duplicate compound_id in index")

    d = artifact_dir(spec)
    d.mkdir(parents = True, exist_ok = True)

    index = pd.DataFrame({
        "compound_id": compound_ids,
        "is_control": is_control if is_control is not None else [False] * len(compound_ids),
    })
    index.to_parquet(d / "index.parquet", index=False)
    np.save(d / "embeddings.npy", matrix)

    payload = {
        "spec": str(spec),
        "n": int(matrix.shape[0]),
        "dim": int(matrix.shape[1]) if matrix.ndim > 1 else 1,
        "dtype": str(matrix.dtype),
        "normalized": False,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        **(meta or {}),
    }
    (d / "meta.json").write_text(json.dumps(payload, indent=2))
    return d


def read_meta(spec: Spec | str) -> dict:
    return json.loads((artifact_dir(spec) / "meta.json").read_text())


def load_embedding(
    spec: Spec | str,
    compound_ids: list[str] | None = None,
    *,
    normalize: bool = False,
    standardize: bool = False,
) -> np.ndarray:
    """
    Load an artifact, reindexed to `compound_ids`.
    """
    d = artifact_dir(spec)
    matrix = np.load(d / "embeddings.npy")
    index = pd.read_parquet(d / "index.parquet")

    if compound_ids is not None:
        pos = pd.Series(np.arange(len(index)), index=index["compound_id"])
        missing = [c for c in compound_ids if c not in pos.index]
        if missing:
            raise KeyError(f"{spec}: {len(missing)} compound_ids absent, e.g. {missing[:5]}")
        matrix = matrix[pos.loc[compound_ids].to_numpy()]

    matrix = matrix.astype(np.float32, copy=True)
    if standardize:
        matrix -= matrix.mean(0, keepdims=True)
        matrix /= matrix.std(0, keepdims=True).clip(min=1e-9)
    if normalize:
        matrix /= np.linalg.norm(matrix, axis=1, keepdims=True).clip(min=1e-9)
    return matrix