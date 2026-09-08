"""
SciPlex-3 (Srivatsan et al. 2020) - scPerturb build SrivatsanTrapnell2020_sciplex3.h5ad
SMILES lookup saved to data/cache/sciplex3_smiles.csv with perturbation, chembl_id, smiles, source
"""

from __future__ import annotations

import time

import anndata as ad
import pandas as pd

from chemembed.config import CACHE, RAW

H5AD = RAW / "sciplex3.h5ad"
SMILES_CACHE = CACHE / "sciplex3_smiles.csv"

DRUG_COL = "perturbation"
CHEMBL_COL = "chembl-ID"
CONTROL_LABELS = {"control", "vehicle", "dmso"}
JUNK_LABELS = {"nan", "none", ""}       

RETRIES = 3
PAUSE = 0.2                             

# manual corrections, just verified last week
MANUAL_SMILES = {
    # PubChem matches by the 2nd word, not the correct  compound
    "INO-1001 (3-Aminobenzamide)": "NC(=O)c1cccc(N)c1",
    "S-Ruxolitinib (INCB018424)": "N#CC[C@@H](C1CCCC1)n1cc(-c2ncnc3[nH]ccc23)cn1",
}

def load_raw(refresh: bool = False) -> pd.DataFrame:
    obs = ad.read_h5ad(H5AD, backed="r").obs
    compounds = _unique_compounds(obs)

    resolved = _resolve(compounds, refresh=refresh)
    compounds["smiles_raw"] = compounds[DRUG_COL].map(resolved)
    compounds.loc[compounds["is_control"], "smiles_raw"] = None

    missing = compounds.loc[~compounds["is_control"] & compounds["smiles_raw"].isna(), DRUG_COL]
    if len(missing):
        print(f"[sciplex3] {len(missing)} unresolved: {missing.tolist()}")
    print(f"[sciplex3] {len(compounds)} rows ({compounds['is_control'].sum()} control)")

    out = compounds.rename(columns={DRUG_COL: "dataset_drug_name"})
    return out[["dataset_drug_name", "smiles_raw", "is_control"]].reset_index(drop=True)


def _unique_compounds(obs: pd.DataFrame) -> pd.DataFrame:
    """1 row per drug name, with its ChEMBL ID"""
    df = pd.DataFrame({
        DRUG_COL: obs[DRUG_COL].astype("string").fillna("").str.strip(),
        "chembl_id": obs[CHEMBL_COL].astype("string").str.strip(),
    })
    df = df[~df[DRUG_COL].str.lower().isin(JUNK_LABELS)]
    df = df.drop_duplicates(DRUG_COL).sort_values(DRUG_COL).reset_index(drop=True)
    df["is_control"] = df[DRUG_COL].str.lower().isin(CONTROL_LABELS)
    return df


def _resolve(compounds: pd.DataFrame, refresh: bool = False) -> dict[str, str]:
    #Return {drug_name: smiles}
    if SMILES_CACHE.exists() and not refresh:
        cache = pd.read_csv(SMILES_CACHE, dtype=str)
    else:
        cache = pd.DataFrame(columns=[DRUG_COL, "chembl_id", "smiles", "source"])
    known = set(cache[DRUG_COL])

    todo = compounds[~compounds["is_control"] & ~compounds[DRUG_COL].isin(known)]
    if len(todo):
        print(f"[sciplex3] resolving {len(todo)} compounds (cached: {len(known)})")
        cache = pd.concat([cache, pd.DataFrame(_lookup_all(todo))], ignore_index=True)
        cache.to_csv(SMILES_CACHE, index=False)
        print(f"[sciplex3] cache -> {SMILES_CACHE}")
        print(cache["source"].value_counts().to_string())

    resolved = cache.dropna(subset=["smiles"])
    return dict(zip(resolved[DRUG_COL], resolved["smiles"]))


def _lookup_all(todo: pd.DataFrame) -> list[dict]:
    from chembl_webresource_client.new_client import new_client
    molecule = new_client.molecule

    rows = []
    for r in todo.itertuples(index=False):
        name = getattr(r, DRUG_COL)
        if name in MANUAL_SMILES:
            rows.append({DRUG_COL: name, "chembl_id": r.chembl_id,
                     "smiles": MANUAL_SMILES[name], "source": "manual"})
            continue
        smiles = _from_chembl(molecule, r.chembl_id)
        source = "chembl"
        if smiles is None:
            smiles = _from_pubchem(name)
            source = "pubchem" if smiles else "none"
        rows.append({DRUG_COL: name, "chembl_id": r.chembl_id, "smiles": smiles, "source": source})
        time.sleep(PAUSE)
    return rows


def _from_chembl(molecule, chembl_id) -> str | None:
    if pd.isna(chembl_id) or not str(chembl_id).strip():
        return None
    for _ in range(RETRIES):
        try:
            hits = list(molecule.filter(molecule_chembl_id=chembl_id).only(["molecule_structures"]))
            if hits and hits[0].get("molecule_structures"):
                return hits[0]["molecule_structures"].get("canonical_smiles")
            return None
        except Exception:
            time.sleep(1)
    return None


def _from_pubchem(name: str) -> str | None:
    import pubchempy as pcp
    for _ in range(RETRIES):
        try:
            hits = pcp.get_compounds(name, "name")
            return hits[0].smiles if hits else None
        except Exception:
            time.sleep(1)
    return None


ANNOTATION_CACHE = CACHE / "sciplex3_annotation.csv"
ANNOTATION_COLS = ["target", "pathway", "pathway_level_1", "pathway_level_2"]

def load_annotation() -> pd.DataFrame:
    """
    SciPlex's own curated drug annotation, one row per drug name.
    chemembed/text/semantic.py for further info
    """
    if ANNOTATION_CACHE.exists():
        return pd.read_csv(ANNOTATION_CACHE, dtype=str)

    obs = ad.read_h5ad(H5AD, backed="r").obs
    df = pd.DataFrame({
        "dataset_drug_name": obs[DRUG_COL].astype("string").fillna("").str.strip(),
        **{c: obs[c].astype("string").str.strip() for c in ANNOTATION_COLS},
    })
    df = df[~df["dataset_drug_name"].str.lower().isin(JUNK_LABELS)]
    df = df.drop_duplicates("dataset_drug_name").reset_index(drop=True)

    df.to_csv(ANNOTATION_CACHE, index=False)
    print(f"[sciplex3] annotation for {len(df)} drugs -> {ANNOTATION_CACHE}")
    return df