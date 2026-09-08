"""
Pipeline: h5ad -> 
    X: embeddings 
    Y: perutbration gene expression responses 

Output: 
    read once per cell line x dose, and cached to data/responses

Y = mean log norm expression - mean of control (same cell line)
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from chemembed.artifacts import load_embedding
from chemembed.compounds import load_all
from chemembed.config import CACHE
from evaluation import config as cfg


# BUILD X AND Y! 

def build_responses(cell_line: str, dose: float, force: bool = False) -> dict:
    """
    Creates the Y matrix per cell line x dose 
    Returns {"Y": (n_compounds, n_genes), "compound_ids": [....], "genes": [....]}
    """
    # formatting stuff 
    path = cfg.RESPONSES / f"{cfg.DATASET}_{cell_line}_{int(dose)}.npz"
    if path.exists() and not force:
        return _read_npz(path)
    import anndata as ad
    import scanpy as sc
    from chemembed.config import RAW
    
    # read in the dataset, filter, log transform, norm, top K genes
    adata = ad.read_h5ad(RAW / f"{cfg.DATASET}.h5ad")
    
    adata.var_names = adata.var["ensembl_id"].astype(str)
    adata.var_names_make_unique()
    
    adata = adata[adata.obs["cell_line"] == cell_line].copy()
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)

    # HVGs from control cells only
    controls = adata.obs["perturbation"].astype(str).str.lower().eq("control")
    hvg = sc.pp.highly_variable_genes(adata[controls].copy(),
                                      n_top_genes = cfg.N_HVG, inplace = False)
    adata = adata[:, hvg["highly_variable"].to_numpy()].copy()

    # control mean 
    obs = adata.obs
    at_time = obs["time"] == cfg.TIME  # TODO: ask - currently picking 24 out of 24 and 72h
    is_control = obs["perturbation"].astype(str).str.lower().eq("control")
    control_mean = np.asarray(adata[is_control & at_time].X.mean(axis = 0)).ravel()
    
    
    # updated expression = treated compounds - control mean (relative to control)
    treated = obs.index[(~is_control) & at_time & (obs["dose_value"] == dose)]
    sub = adata[treated]
    names = sub.obs["perturbation"].astype(str).str.strip()
    rows, kept_names = [], []
    for name, idx in names.groupby(names).groups.items():
        if len(idx) < cfg.MIN_CELLS_PER_CONDITION:
            continue
        mean = np.asarray(sub[idx].X.mean(axis = 0)).ravel()
        rows.append(mean - control_mean)
        kept_names.append(name)

    Y = np.vstack(rows).astype(np.float32) # Y.shape: #compounds x #genes
    ids = _names_to_compound_ids(kept_names)

    keep = [i for i, c in enumerate(ids) if c is not None]
    Y, ids = Y[keep], [ids[i] for i in keep]

    # edge case, multiple compounds map to the same structure 
    if len(set(ids)) != len(ids):
        order = pd.unique(pd.Series(ids))
        Y = np.vstack([Y[[i for i, c in enumerate(ids) if c == cid]].mean(axis = 0)
                       for cid in order]).astype(np.float32)
        print(f"[responses] collapsed {len(ids)} rows -> {len(order)} unique compounds")
        ids = list(order)
    
    np.savez_compressed(path, Y = Y,
                        compound_ids = np.array(ids),
                        genes = adata.var_names.to_numpy(dtype = str),
                        control_mean = control_mean.astype(np.float32))
    print(f"[responses] {cell_line} @ {dose:g}nM: {Y.shape} -> {path}")
    return _read_npz(path)


def _read_npz(path) -> dict:
    """helper func"""
    z = np.load(path, allow_pickle = False)
    out = {"Y": z["Y"], "compound_ids": list(z["compound_ids"]), "genes": list(z["genes"])}
    # older caches predate control_mean; rebuild them to get absolute metrics
    out["control_mean"] = z["control_mean"] if "control_mean" in z.files else None
    return out


def _names_to_compound_ids(names: list[str]) -> list[str | None]:
    """ maps obs drug name -> InChIKey, done with the dataset specific compound table"""
    from chemembed.config import COMPOUNDS
    table = pd.read_csv(COMPOUNDS / f"{cfg.DATASET}.csv")
    lookup = dict(zip(table["dataset_drug_name"], table["compound_id"]))
    missing = [n for n in names if n not in lookup]
    if missing:
        print(f"[responses] {len(missing)} names not in the compound table, dropped: "
              f"{missing[:5]}")
    mapping = [lookup.get(n) for n in names]
    return mapping

def load_xy(spec: str, cell_line: str, dose: float = cfg.PRIMARY_DOSE):
    """ returns (X, Y, compound_ids, control_mean), rows aligned by compound_id """
    resp = build_responses(cell_line, dose)
    ids = [c for c in resp["compound_ids"] if not c.startswith(cfg.CONTROL_ID_PREFIX)]
    wanted = set(ids)
    keep = [i for i, c in enumerate(resp["compound_ids"]) if c in wanted]
    Y = resp["Y"][keep]
    X = load_embedding(spec, ids)
    return X, Y, ids, resp["control_mean"]

def load_xy_doses(spec: str, cell_line: str, doses: tuple = cfg.DOSES,
                  interaction: bool = True):
    """
    All doses stacked, with dose as a feature. Adapted from CIPHER's
    add_dose_feature. X = [E, log_dose, E * log_dose]
    """
    Xs, Ys, ids_all, dose_col, control_means = [], [], [], [], []

    for dose in doses:
        resp = build_responses(cell_line, dose)
        ids = [c for c in resp["compound_ids"]
               if not c.startswith(cfg.CONTROL_ID_PREFIX)]
        wanted = set(ids)
        keep = [i for i, c in enumerate(resp["compound_ids"]) if c in wanted]

        Ys.append(resp["Y"][keep])
        Xs.append(load_embedding(spec, ids))
        ids_all += ids
        dose_col += [dose] * len(ids)
        control_means.append(resp["control_mean"])

    E = np.vstack(Xs)
    Y = np.vstack(Ys).astype(np.float32)

    # log dose, standardised. Dose is a design variable known for test rows too,
    # so standardising over all rows is not outcome leakage.
    log_dose = np.log10(np.asarray(dose_col, dtype = np.float64) + 1.0)[:, None]
    log_dose = ((log_dose - log_dose.mean()) / (log_dose.std() + 1e-8)).astype(E.dtype)

    parts = [E, log_dose] + ([E * log_dose] if interaction else [])
    X = np.hstack(parts).astype(np.float32)

    return X, Y, ids_all, control_means[0]


# SPLITS 
def make_folds(ids: list[str], kind: str, n_folds: int = cfg.N_FOLDS):
    """
    returns a list of (train_idx, test_idx) for each fold, compound-level.

    When ids repeat - which happens once doses are stacked - splitting on row
    index would put the same compound in train and test. Grouping prevents that.
    """
    from sklearn.model_selection import GroupKFold, KFold
    seed = 17
    rng = np.random.default_rng(seed)

    if kind == "random":
        if len(set(ids)) == len(ids):
            kf = KFold(n_splits = n_folds, shuffle = True, random_state = seed)
            return list(kf.split(np.arange(len(ids))))
        groups = np.asarray(ids)                      # duplicates -> group by compound

    elif kind == "scaffold":
        groups = scaffold_groups(ids)

    else:
        raise ValueError(f"split kind not known {kind!r}")

    # GroupKFold is deterministic; shuffle the group order so folds vary with seed
    order = {g: i for i, g in enumerate(rng.permutation(pd.unique(groups)))}
    shuffled = np.array([order[g] for g in groups])
    return list(GroupKFold(n_splits = n_folds).split(np.arange(len(ids)),
                                                     groups = shuffled))
    

def scaffold_groups(ids: list[str]) -> np.ndarray:
    """ looks up murcko scaffolds """
    descriptions = json.loads((CACHE / "descriptions.json").read_text())
    scaffolds = []
    for c in ids:
        text = descriptions.get(c, {}).get("murcko_scaffold", "")
        scaffolds.append(text.replace("Murcko scaffold:", "").strip() or c)
    return np.array(scaffolds)



# EVAL METRICS 
def evaluate(Y_true, Y_pred, control_mean = None, top_n: int = cfg.TOP_N_DEGS) -> dict:
    """
    metrics are calculated per compound, then averaged over compounds.
    
    Y: delta expression, shape: (n_compounds, n_genes)
    DE genes: top delta genes of the OBSERVED response, per compound.
    """

    n = len(Y_true)
    keys = ("pearson", "pearson_top", "r2", "mse", "mae", "direction", "overlap_top")
    if control_mean is not None:
        keys = keys + ("r2_abs", "r2_abs_top")
    per = {k: np.zeros(n) for k in keys}

    for i in range(n):
        true, pred = Y_true[i], Y_pred[i]
        top = np.argsort(-np.abs(true))[:top_n]

        per["pearson"][i] = pearson(true, pred)
        per["pearson_top"][i] = pearson(true[top], pred[top])
        per["r2"][i] = r2(true, pred)
        per["mse"][i] = np.mean((true - pred) ** 2)
        per["mae"][i] = np.mean(np.abs(true - pred))
        per["direction"][i] = np.mean(np.sign(true) == np.sign(pred))
        pred_top = np.argsort(-np.abs(pred))[:top_n]
        per["overlap_top"][i] = len(np.intersect1d(top, pred_top)) / top_n

        if control_mean is not None:
            abs_true, abs_pred = true + control_mean, pred + control_mean
            per["r2_abs"][i] = r2(abs_true, abs_pred)
            per["r2_abs_top"][i] = r2(abs_true[top], abs_pred[top])

    out = {k: float(np.mean(v)) for k, v in per.items()}
    out["r2_gene_sd"] = float(r2(Y_true.std(axis = 0), Y_pred.std(axis = 0)))
    return out


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a - a.mean(), b - b.mean()
    denom = np.sqrt((a * a).sum() * (b * b).sum())
    return float(a @ b / denom) if denom > 1e-17 else 0.0


def r2(true: np.ndarray, pred: np.ndarray) -> float:
    ss_res = np.sum((true - pred) ** 2)
    ss_total = np.sum((true - true.mean()) ** 2)
    return float(1 - ss_res / ss_total) if ss_total > 1e-17 else 0.0