
# Shared backbone across generator scripts.

from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
from chemembed.artifacts import write_artifact
from chemembed.compounds import load_all
from chemembed.registry import Spec


def parse_spec() -> Spec:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    return Spec.parse(ap.parse_args().spec)


def load_inputs() -> pd.DataFrame:
    df = load_all()
    keep = ~df["is_control"].astype(bool) & df["smiles_canonical"].notna()
    return df[keep].reset_index(drop=True)


def finish(spec: Spec, compound_ids: list[str], matrix: np.ndarray, meta: dict | None = None):
    all_compounds = load_all()
    controls = all_compounds.loc[all_compounds["is_control"].astype(bool), "compound_id"].tolist()
    if controls:
        matrix = np.vstack([matrix, np.zeros((len(controls), matrix.shape[1]), matrix.dtype)])
        compound_ids = list(compound_ids) + controls
    path = write_artifact(
        spec, compound_ids, matrix,
        is_control=[False] * (len(compound_ids) - len(controls)) + [True] * len(controls),
        meta=meta,
    )
    print(f"[done] {spec}: {matrix.shape} -> {path}")