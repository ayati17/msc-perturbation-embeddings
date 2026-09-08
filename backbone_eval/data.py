"""
Data Preprocessing for Cellflow 
h5ad -> the two AnnData objects CellFlow needs, plus compound-level folds.

CellFlow trains on control + seen compounds and predicts held-out compounds from
control cells, so the split is by compound and controls appear in BOTH objects.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from backbone_eval import config as cfg


def load_adata():
    import scanpy as sc
    from chemembed.config import COMPOUNDS, RAW

    cache = cfg.RESULTS / "sciplex3_preprocessed.h5ad"

    if cache.exists():
        print(f"[data] loading cached {cache}", flush=True)
        return sc.read_h5ad(cache)

    print("[data] preprocessing from raw...", flush=True)

    adata = sc.read_h5ad(RAW / f"{cfg.DATASET}.h5ad")
    
    adata.var_names = adata.var["ensembl_id"].astype(str)
    adata.var_names_make_unique()
    
    adata.obs["perturbation"] = adata.obs["perturbation"].astype(str).str.strip()

    table = pd.read_csv(COMPOUNDS / f"{cfg.DATASET}.csv")
    known = set(table["dataset_drug_name"].str.strip())
    is_control = adata.obs["pathway_level_1"].eq("Vehicle")

    keep = (
        adata.obs["time"].eq(cfg.TIME)
        & adata.obs["cell_line"].notna()
        & (adata.obs["perturbation"].isin(known) | is_control)
    )
    adata = adata[keep.values].copy()

    adata.obs["is_control"] = adata.obs["pathway_level_1"].eq("Vehicle").values
    adata.obs["log_dose"] = np.log10(
        adata.obs["dose_value"].clip(lower=1.0)
    )
    adata.obs.loc[adata.obs["is_control"], "log_dose"] = 0.0

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=cfg.N_HVG,
        flavor="seurat",
        batch_key="cell_line",
    )

    adata = adata[:, adata.var["highly_variable"]].copy()

    print(f"[data] saving cache -> {cache}", flush=True)
    adata.write_h5ad(cache)

    return adata


def compound_folds(adata, n_folds: int = cfg.N_FOLDS, seed: int = cfg.SEED):
    """ {fold: set of held-out compound names}, never control"""
    from sklearn.model_selection import KFold

    drugs = sorted(set(adata.obs.loc[~adata.obs["is_control"], "perturbation"]))
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    return {f: {drugs[i] for i in test} for f, (_, test) in enumerate(kf.split(drugs))}


def split(adata, test_drugs: set):
    """
    (train, test) AnnData with a shared PCA basis fitted on train compounds. 
    """
    from cellflow.preprocessing import centered_pca, project_pca

    held_out = adata.obs["perturbation"].isin(test_drugs) & ~adata.obs["is_control"]
    train = adata[~held_out.values].copy()
    test = adata[(held_out | adata.obs["is_control"]).values].copy()

    centered_pca(train, n_comps=cfg.N_PCA, keep_centered_data=False)
    project_pca(test, ref_adata=train)

    print(f"[split] train {train.n_obs:,} cells | test {test.n_obs:,} cells "
          f"({len(test_drugs)} held-out compounds)")
    return train, test