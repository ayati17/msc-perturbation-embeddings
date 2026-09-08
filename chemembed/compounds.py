"""
Dataset -> compounds table.
"""
from __future__ import annotations
import pandas as pd
from chemembed.config import COMPOUNDS

COLUMNS = [
    "compound_id",
    "dataset",
    "dataset_drug_name",
    "smiles_raw",
    "smiles_canonical",
    "is_control",
    "pubchem_cid",
]


def build_compounds(dataset: str) -> pd.DataFrame:
    # Imported here, not at module scope: generator environments import load_all()
    # from this file and must not be made to install rdkit and anndata for it.
    from chemembed.canonicalize import standardize
    from chemembed.datasets import LOADERS

    if dataset not in LOADERS:
        raise ValueError(f"unknown dataset {dataset!r}, have {list(LOADERS)}")

    raw = LOADERS[dataset]()
    rows, failed = [], []
    for r in raw.itertuples(index=False):
        if r.is_control:
            rows.append({
                "compound_id": f"CONTROL_{dataset.upper()}",
                "smiles_canonical": None,
            } | _base(dataset, r))
            continue
        if pd.isna(r.smiles_raw):
            failed.append(r.dataset_drug_name)
            continue
        std = standardize(r.smiles_raw)
        if std is None:
            failed.append(r.dataset_drug_name)
            continue
        rows.append({
            "compound_id": std.inchikey,
            "smiles_canonical": std.smiles,
        } | _base(dataset, r))

    df = pd.DataFrame(rows).reindex(columns=COLUMNS)
    df = df.drop_duplicates(["compound_id", "dataset_drug_name"]).reset_index(drop=True)
    df.to_csv(COMPOUNDS / f"{dataset}.csv", index=False)

    if failed:
        print(f"[compounds] {dataset}: {len(failed)} without a usable structure, "
              f"e.g. {failed[:5]}")
    collisions = len(df) - df["compound_id"].nunique()
    if collisions:
        print(f"[compounds] {dataset}: {collisions} name(s) share a compound_id with another")
    print(f"[compounds] {dataset}: {len(df)} names, {df['compound_id'].nunique()} structures "
          f"-> {COMPOUNDS / f'{dataset}.csv'}")

    _rebuild_all()
    return df


def _base(dataset: str, r) -> dict:
    return {
        "dataset": dataset,
        "dataset_drug_name": r.dataset_drug_name,
        "smiles_raw": r.smiles_raw,
        "is_control": bool(r.is_control),
        "pubchem_cid": pd.NA,
    }


def _rebuild_all() -> pd.DataFrame:
    parts = [pd.read_csv(p) for p in sorted(COMPOUNDS.glob("*.csv")) if p.name != "all.csv"]
    allc = pd.concat(parts).drop_duplicates("compound_id").reset_index(drop=True)
    allc.to_csv(COMPOUNDS / "all.csv", index=False)
    print(f"[compounds] all.csv: {len(allc)} structures across {allc['dataset'].nunique()} dataset(s)")
    return allc


def duplicate_names(dataset: str) -> pd.DataFrame:
    """Drug names that collapsed onto the same structure - worth eyeballing once."""
    df = pd.read_csv(COMPOUNDS / f"{dataset}.csv")
    dupes = df[df.duplicated("compound_id", keep=False)]
    return dupes.sort_values("compound_id")[["compound_id", "dataset_drug_name", "smiles_canonical"]]


def load_all() -> pd.DataFrame:
    return pd.read_csv(COMPOUNDS / "all.csv")