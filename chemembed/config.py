import os
from pathlib import Path

ROOT = Path(os.environ.get("CHEMEMBED_ROOT", Path(__file__).resolve().parents[1]))

DATA = ROOT / "data"
RAW = DATA / "raw"            # h5ad files, just read-only # TODO: add in the datasets here 
COMPOUNDS = DATA / "compounds" # <dataset>.csv, all.csv
CACHE = DATA / "cache" # pubchem lookups, rendered description text
ARTIFACTS = DATA / "artifacts" # <spec>/{embeddings.npy,index.parquet,meta.json}

# 1 virtualenv per environment name: <ENVS>/<name>/bin/python
ENVS = Path(os.environ.get("CHEMEMBED_ENVS", ROOT / ".envs"))

for _p in (RAW, COMPOUNDS, CACHE, ARTIFACTS):
    _p.mkdir(parents = True, exist_ok = True)