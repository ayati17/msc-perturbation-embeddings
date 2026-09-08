
"""
PubChem lookups: CID and IUPAC name
env: core

Output: data cached to data/cache/pubchem.json, keyed by compound_id.
Steps:
1. Run once after build_compounds(), BEFORE build_descriptions() and any text spec:
    build_pubchem()             
    
Other info: 
PubChem
 - 5 requests/sec allowed in pubchem 
 - IUPAC name is one of six structural renderings in text/descriptions.py, and PubChem covers it completely (187/187 on SciPlex-3).
 - originally tried PubChem's compound description, but returned 0% DrugBank, 0% LiverTox and 79% ChEBI (just structural restatements)

updated semantic text - ChEMBL mechanism of action + SciPlex's annotations inchemembed/text/semantic.py.
"""

from __future__ import annotations

import json
import time

import pandas as pd
import requests

from chemembed.compounds import load_all
from chemembed.config import CACHE

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CACHE_FILE = CACHE / "pubchem.json"

MIN_INTERVAL = 0.25                                  

_session = requests.Session()
_last_call = 0.0


def build_pubchem(limit: int | None = None, force: bool = False) -> dict:
    """Fill the cache for every compound in all.csv. Returns {compound_id: record}."""
    cache = {} if force else _read_cache()

    compounds = load_all()
    todo = compounds[~compounds["is_control"].astype(bool)]
    todo = todo[~todo["compound_id"].isin(cache)]
    if limit:
        todo = todo.head(limit)

    print(f"[pubchem] {len(todo)} to fetch, {len(cache)} cached")
    for i, r in enumerate(todo.itertuples(index=False), 1):
        cache[r.compound_id] = _fetch(r.smiles_canonical, r.dataset_drug_name)
        if i % 25 == 0:
            print(f"[pubchem] {i}/{len(todo)}")
            _write_cache(cache)

    _write_cache(cache)
    _report(cache)
    return cache


def load_iupac_names() -> dict[str, str]:
    return {k: (v.get("iupac_name") or "") for k, v in _read_cache().items()}


def load_cids() -> dict[str, int | None]:
    return {k: v.get("cid") for k, v in _read_cache().items()}


# for 1 chemical compound 

def _fetch(smiles: str, name: str) -> dict:
    """SMILES identifies the molecule; name is only a fallback for CID lookup."""
    cid = _cid_from_smiles(smiles) or _cid_from_name(name)
    if cid is None:
        return {"cid": None, "iupac_name": None}
    return {"cid": cid, "iupac_name": _iupac_name(cid)}


def _cid_from_smiles(smiles: str) -> int | None:
    data = _request("POST", f"{PUG}/compound/smiles/cids/JSON", data={"smiles": smiles})
    return _first_cid(data)


def _cid_from_name(name: str) -> int | None:
    data = _request("GET", f"{PUG}/compound/name/{requests.utils.quote(name)}/cids/JSON")
    return _first_cid(data)


def _iupac_name(cid: int) -> str | None:
    data = _request("GET", f"{PUG}/compound/cid/{cid}/property/IUPACName/JSON")
    props = (data or {}).get("PropertyTable", {}).get("Properties", [])
    return props[0].get("IUPACName") if props else None

def _first_cid(data: dict | None) -> int | None:
    cids = (data or {}).get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def _request(method: str, url: str, **kwargs) -> dict | None:
    global _last_call
    time.sleep(max(0.0, MIN_INTERVAL - (time.time() - _last_call)))
    _last_call = time.time()
    try:
        response = _session.request(method, url, timeout=30, **kwargs)
        return response.json() if response.ok else None
    except Exception as e:
        print(f"[pubchem] {url} failed: {e}")
        return None


def _read_cache() -> dict:
    return json.loads(CACHE_FILE.read_text()) if CACHE_FILE.exists() else {}


def _write_cache(cache: dict) -> None:
    CACHE_FILE.write_text(json.dumps(cache, indent=2))


def _report(cache: dict) -> None:
    missing_cid = sum(1 for v in cache.values() if not v.get("cid"))
    missing_iupac = sum(1 for v in cache.values() if not v.get("iupac_name"))
    print(f"[pubchem] {len(cache)} compounds -> {CACHE_FILE}")
    print(f"[pubchem] missing CID: {missing_cid}, missing IUPAC name: {missing_iupac}")