"""
Semantic (non-structural) text for molecular function 
env: core
notes for further analysis - how much did dataset-specific curation add over public biological data

chEMBL mechanism where available + SciPlex's own target pathway annotation.
Output: 
    data/cache/semantic.json    {compound_id: {chembl, curated, source, has_mechanism}}

Steps: 
1. build_semantic(), to be run after build_compounds():
"""

from __future__ import annotations

import json
import time

import pandas as pd


from chemembed.config import CACHE, COMPOUNDS

CACHE_FILE = CACHE / "semantic.json"
MECHANISM_CACHE = CACHE / "chembl_mechanisms.json"
MECHANISM_URL = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
PAUSE = 0.2


def build_semantic(dataset: str = "sciplex3") -> dict:
    compounds = pd.read_csv(COMPOUNDS / f"{dataset}.csv")
    compounds = compounds[~compounds["is_control"].astype(bool)]

    chembl_ids = _chembl_ids(dataset)
    mechanisms = _mechanisms(sorted(set(chembl_ids.values()) - {None}))
    annotation = _annotation(dataset)

    out = {}
    for r in compounds.itertuples(index=False):
        name = r.dataset_drug_name
        found = mechanisms.get(chembl_ids.get(name), [])
        note = annotation.get(name, {})

        chembl_text = _render_mechanism(name, found)
        curated_text = " ".join(t for t in (chembl_text or f"{name}.",
                                            _render_annotation(note)) if t).strip()

        out[r.compound_id] = {
            "chembl": chembl_text,
            "curated": curated_text,
            "source": "chembl" if found else ("annotation" if note else "none"),
            "has_mechanism": bool(found),
        }

    CACHE_FILE.write_text(json.dumps(out, indent=2))
    _report(out)
    return out


def load_semantic(variant: str) -> dict[str, str]:
    """variant: 'chembl' | 'curated'"""
    return {k: v.get(variant, "") for k, v in json.loads(CACHE_FILE.read_text()).items()}

def load_sources() -> dict[str, str]:
    return {k: v.get("source", "none")
            for k, v in json.loads(CACHE_FILE.read_text()).items()}
    
def covered_compounds() -> list[str]:
    """compound_ids with a real ChEMBL mechanism - the subset for the clean comparison."""
    cache = json.loads(CACHE_FILE.read_text())
    return [k for k, v in cache.items() if v.get("has_mechanism")]


def _render_mechanism(name: str, mechanisms: list[str]) -> str:
    if not mechanisms:
        return ""
    return f"{name}. Mechanism of action: {'; '.join(mechanisms)}."


def _render_annotation(note: dict) -> str:
    parts = []
    if note.get("target"):
        parts.append(f"Annotated target: {note['target']}.")
    pathway = note.get("pathway_level_1") or note.get("pathway")
    if pathway:
        detail = note.get("pathway_level_2")
        parts.append(f"Pathway: {pathway}{f' ({detail})' if detail else ''}.")
    return " ".join(parts)


def _chembl_ids(dataset: str) -> dict[str, str | None]:
    """drug name -> first ChEMBL ID, from the SMILES resolution cache."""
    smiles = pd.read_csv(CACHE / f"{dataset}_smiles.csv", dtype=str)
    ids = smiles["chembl_id"].str.split(";").str[0]
    return dict(zip(smiles["perturbation"], ids.where(ids.notna(), None)))


def _annotation(dataset: str) -> dict[str, dict]:
    path = CACHE / f"{dataset}_annotation.csv"
    if not path.exists():
        print(f"[semantic] no annotation at {path}; run load_annotation() first")
        return {}
    df = pd.read_csv(path, dtype=str).fillna("")
    return {r.pop("dataset_drug_name"): r for r in df.to_dict("records")}


def _mechanisms(chembl_ids: list[str]) -> dict[str, list[str]]:
    """ ChEMBL MOA data, cached - one request per compound."""
    cache = json.loads(MECHANISM_CACHE.read_text()) if MECHANISM_CACHE.exists() else {}
    todo = [c for c in chembl_ids if c not in cache]
    print(f"[semantic] {len(todo)} ChEMBL lookups ({len(cache)} cached)")

    for i, chembl_id in enumerate(todo, 1):
        data = _get(MECHANISM_URL, params={"molecule_chembl_id": chembl_id})
        cache[chembl_id] = [m["mechanism_of_action"]
                            for m in (data or {}).get("mechanisms", [])
                            if m.get("mechanism_of_action")]
        if i % 25 == 0:
            MECHANISM_CACHE.write_text(json.dumps(cache, indent=2))

    MECHANISM_CACHE.write_text(json.dumps(cache, indent=2))
    return cache


def _get(url: str, **kwargs) -> dict | None:
    import requests
    time.sleep(PAUSE)
    try:
        response = requests.get(url, timeout=30, **kwargs)
        return response.json() if response.ok else None
    except Exception as e:
        print(f"[semantic] {url} failed: {e}")
        return None


def _report(out: dict) -> None:
    counts = pd.Series([v["source"] for v in out.values()]).value_counts()
    empty_chembl = sum(1 for v in out.values() if not v["chembl"])
    empty_curated = sum(1 for v in out.values() if not v["curated"].strip(". "))
    print(f"[semantic] {len(out)} compounds -> {CACHE_FILE}")
    print(f"[semantic] source:\n{counts.to_string()}")
    print(f"[semantic] empty text - chembl: {empty_chembl}, curated: {empty_curated}")