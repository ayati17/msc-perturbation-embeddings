from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple
import numpy as np
import pandas as pd


@dataclass
class Config:
    # data and config setup 
    adata_path: str = "data/raw/sciplex3.h5ad"
    work_dir: str = "data/results/cipher"
    dataset: str = "sciplex3"

    n_hvg: int = 2000
    min_cells_per_group: int = 5
    min_control_cells: int = 50

    residual_formulation: bool = True
    dose_interaction: bool = True
    embeddings: Tuple[str, ...] = () # list of embeddings to run (empty means all files) 
    
    random_state: int = 17
    n_outer_folds: int = 5
    n_inner_folds: int = 4

    # hyerpaparameter tuning! 
    # 1. shrinkage - controls how strongly sigma is regularised. larger the val, more regularised sigma (less specific)
    # 2. lamba - ridge regression lambda 
    # 3. gamma - the exponent, controls how penalty would be distributed across the covariance modes
    shrinkage_grid: Tuple[float, ...] = (0.01, 0.05, 0.15, 0.4)
    lambda_grid: Tuple[float, ...] = tuple(np.logspace(-1, 7, 25))
    gamma_grid: Tuple[float, ...] = (-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)

    def __post_init__(self):
        Path(self.work_dir).mkdir(parents = True, exist_ok = True)


# =========
# SECTION 0 
# 1. Constructs 
#       Delta X = <x>_u - <x>_0 
#       for every (compound, dose, cell line) group combination 
# 2. Keep track of control cells for further calculation
# 3. Group folds for unseen compound holdout
# ==========

def build_response_matrix(adata, cfg: Config):
    """
    Inputs: 
    adata: sciplex3 dataset 
    
    Outputs: 
    (DX, META, X, obs, gene_names, is_control).
    DX[r]  = mean expression of combo group r  -  mean expression of its control
    META   = one row per group: perturbation, dose_value, cell_line, n_cells
    """
    import scanpy as sc

    # Normalizes cells, keeps top 2000 
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = cfg.n_hvg, flavor = "seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Calculate all the control cells per cell line 
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X = X.astype(np.float32)
    obs = adata.obs
    is_control = (obs["perturbation"] == "control").values
    # TODO: clean this up 
    # ========
    gene_names = np.asarray(adata.var_names)
    for cand in ("gene_short_name", "gene_symbol", "gene_name", "symbol", "feature_name"): 
        if cand in adata.var.columns and adata.var[cand].notna().all():
            gene_names = np.asarray(adata.var[cand])
            break
    # ========
    
    ctrl_counts = obs.loc[is_control, "cell_line"].value_counts()
    cell_lines = ctrl_counts[ctrl_counts >= cfg.min_control_cells].index.tolist()
    print("cell lines kept:", cell_lines)
    print("control cells per line:", {cl: int(ctrl_counts[cl]) for cl in cell_lines})

    # Calculate the control mean per cell line 
    control_mean = {cl: X[is_control & (obs["cell_line"] == cl).values].mean(0) for cl in cell_lines}

    # Perutbration cells (everything that's not control)
    keep = (~is_control) & obs["cell_line"].isin(cell_lines).values
    frame = obs.loc[keep, ["perturbation", "dose_value", "cell_line"]].copy()
    frame["row"] = np.where(keep)[0]

    # Creation of groups 
    rows, meta = [], []
    for (pert, dose, cl), sub in frame.groupby(["perturbation", "dose_value", "cell_line"], observed = True):
        idx = sub["row"].values
        if len(idx) < cfg.min_cells_per_group:
            continue
        rows.append(X[idx].mean(0) - control_mean[cl])
        meta.append(dict(perturbation = str(pert).strip(), 
                         dose_value = float(dose), 
                         cell_line = str(cl), 
                         n_cells = int(len(idx))))

    DX = np.vstack(rows).astype(np.float32)    # Matrix Dims: num perturbation groups x num genes (2000)
    META = pd.DataFrame(meta) # experiment details per row of DX 
    print(f"pseudobulk responses: {DX.shape}, median cells/group = " f"{META['n_cells'].median():.0f}") # TODO: clean this up 
    return DX, META, X, obs, gene_names, is_control

def assign_outer_folds(META: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    # Setup folds by group so that unseen compounds are held out entirely into either training or test datasets

    from sklearn.model_selection import GroupKFold

    META = META.copy()
    META["fold"] = -1
    group_kf = GroupKFold(n_splits = cfg.n_outer_folds)
    # Assign folds 
    for f, (_, te) in enumerate(group_kf.split(META, groups = META["perturbation"])):
        META.loc[META.index[te], "fold"] = f
    return META


# =========
# SECTION 1
# eigendecomp of the control covariance for a given cell line 
# 1. calculates the sample estimate of the covariance matrix 
# 2. calculates & stabilizes + sorts the eigenvalues and eigenvectors of the covariance matrix
# 3. effective dimensionality + analysis of how variance is concentrated in the top modes 
# =========
class FluctuationBasis:
    
    def __init__(self, Xc: np.ndarray):
        # control cells - centering 
        Xc = Xc - Xc.mean(0, keepdims = True)
        self.n_cells, self.G = Xc.shape
        
        # S ~= the sigma matrix, calc eigenvalues and eigenvectors
        # sort them into descending order 
        pairwise = (Xc.T @ Xc)
        S =  pairwise / max(self.n_cells - 1, 1)
        eigenvals, eigenvecs = np.linalg.eigh(S.astype(np.float64))
        self.eigenvals = eigenvals[::-1].copy()          
        self.eigenvecs = eigenvecs[:, ::-1].copy()       
        
        # calc average eigenvalue
        trace = self.eigenvals.sum()
        self.trace_over_G = float(trace / self.G)

    def shrink_eigenvals(self, shrinkage: float) -> np.ndarray:
        # stabalise eigenvalues
        lambda_is = ((1.0 - shrinkage) * self.eigenvals) + (shrinkage * self.trace_over_G)
        just_in_case = np.maximum(lambda_is, 1e-12)
        return just_in_case

    def to_modes(self, A: np.ndarray) -> np.ndarray:
        # map rows of A (n x G) to the eigenbasis of covariance matrix 
        mode_space = A @ self.eigenvecs
        return mode_space

    def from_modes(self, C: np.ndarray) -> np.ndarray:
        # map mode coefficients back to gene space 
        gene_space = C @ self.eigenvecs.T
        return gene_space

    def effective_rank(self, shrinkage: float = 0.0) -> float:
        # check how much variance can be attributed to each specific mode
        lambda_is = self.shrink_eigenvals(shrinkage)
        p_is = lambda_is / lambda_is.sum()
        
        # measure of effective dimensionality 
        expo_entropy = (-(p_is * np.log(p_is + 1e-300)).sum())
        return float(np.exp(expo_entropy))

    def variance_dict(self, ks: Sequence[int] = (1, 3, 10, 30, 100, 200, 500)) -> Dict[int, float]:
        total_var = self.eigenvals.sum()
        var_dict = {k: float(self.eigenvals[:k].sum() / total_var) for k in ks}
        return var_dict

def build_bases(X, obs, is_control, cell_lines) -> Dict[str, FluctuationBasis]:
    bases = {}
    for cl in cell_lines:
        mask = is_control & (obs["cell_line"] == cl).values
        bases[cl] = FluctuationBasis(X[mask])
        b = bases[cl]
        print(f"{cl}: n_ctrl={b.n_cells}, effective_rank = {b.effective_rank(0.05):.1f}, "
              f"variance in the top modes {b.variance_dict()}")
    return bases


# =========
# SECTION 2
# mode specific ridge penalty 
# 1. calculates for each i alpha_i = lambda * (lambda_bar / lambda_i) ** 2xgamma
# 2. different values of gamma for hyperparam optimisation 
# =========
def mode_penalties(spectrum: np.ndarray, lam: float, gamma: float) -> np.ndarray:
    avg_variance = float(spectrum.mean())
    penalties = lam * (avg_variance / spectrum) ** (2.0 * gamma)
    return penalties


# =========
# SECTION 3
# fits embedding to response relationship w mode specific ridge penaltie s w 1 svd
# 1. perform SVD decomp of the training embeddings matrix, reuse that for all modes and hyerparams 
# 2. make predictions 
# 3. calc the forces 
# ========

class ModeRidge:
    def __init__(self, basis: FluctuationBasis, shrinkage: float):
        self.basis = basis
        self.spectrum = basis.shrink_eigenvals(shrinkage)

    def fit(self, E_tr: np.ndarray, DX_tr: np.ndarray):
        # standardise the chem embeddings 
        self.mu_E = E_tr.mean(0, keepdims = True)
        self.sd_E = E_tr.std(0, keepdims = True) + 1e-8
        Ec = (E_tr - self.mu_E) / self.sd_E

        # setup Y, standatidise 
        Y = self.basis.to_modes(DX_tr)            
        self.mu_Y = Y.mean(0, keepdims = True)
        Yc = Y - self.mu_Y
        
        # SVD of chem. embeddings 
        U, s, Vt = np.linalg.svd(Ec, full_matrices = False)
        keep = s > s.max() * 1e-10
        self.U, self.s, self.Vt = U[:, keep], s[keep], Vt[keep]
        self.UtY = self.U.T @ Yc                  
        return self

    def predict(self, E_te: np.ndarray, lam: float, gamma: float) -> np.ndarray:
        Ec = (E_te - self.mu_E) / self.sd_E
        P = Ec @ self.Vt.T                        
        alpha = mode_penalties(self.spectrum, lam, gamma)          
        F = self.s[:, None] / (self.s[:, None] ** 2 + alpha[None, :])  
        Y_hat = P @ (F * self.UtY) + self.mu_Y    # n_te x G  (mode space)
        Y_hat_converted = self.basis.from_modes(Y_hat) # n_te x G  (gene space)
        return Y_hat_converted  

    def predict_grid(self, E_te, 
                     lambdas: Iterable[float], gammas: Iterable[float],
                     in_modes: bool = False):
        # repeats predict for all lambda, gamma combinations using one SVD 
        Ec = (E_te - self.mu_E) / self.sd_E
        P = Ec @ self.Vt.T
        for gamma in gammas:
            for lam in lambdas:
                alpha = mode_penalties(self.spectrum, lam, gamma)
                F = self.s[:, None] / (self.s[:, None] ** 2 + alpha[None, :])
                Y_hat = P @ (F * self.UtY) + self.mu_Y
                yield lam, gamma, (Y_hat if in_modes else self.basis.from_modes(Y_hat))

    def implied_force(self, E_te: np.ndarray, lam: float, gamma: float) -> np.ndarray:
        # given sigma, calculates model's implied force u_hat 
        """u_hat itself, for interpretation: which genes is the drug forcing?"""
        Ec = (E_te - self.mu_E) / self.sd_E
        P = Ec @ self.Vt.T
        alpha = mode_penalties(self.spectrum, lam, gamma)
        F = self.s[:, None] / (self.s[:, None] ** 2 + alpha[None, :])
        Y_hat = P @ (F * self.UtY) + self.mu_Y
        u_hat = self.basis.from_modes(Y_hat / self.spectrum[None, :])
        return u_hat


def observed_force(basis: FluctuationBasis, DX: np.ndarray, shrinkage: float) -> np.ndarray:
    # calculates real data's force u_obs = sigma^-1 * delta X
    lam = basis.shrink_eigenvals(shrinkage)
    u_obs = basis.from_modes(basis.to_modes(DX) / lam[None, :])
    return u_obs


# =========
# SECTION 4
# 1. setup baseline for comparison = avg response for a given cell line + factor in the dosage 
# 2 chemical embedding permutation, testing if the embedings contain any useful info 
# =========+=
class StratumMean:
    def __init__(self, keys: Sequence[str] = ("cell_line", "dose_value")):
        self.keys = list(keys)

    def fit(self, DX: np.ndarray, META: pd.DataFrame, train_mask: np.ndarray):
        sub = META.loc[train_mask]
        self.table = {k: DX[np.asarray(g.index)].mean(0)
                      for k, g in sub.groupby(self.keys, observed = True)}
        self.global_mean = DX[train_mask].mean(0)
        return self

    def predict(self, META: pd.DataFrame, idx: np.ndarray) -> np.ndarray:
        out = np.zeros((len(idx), self.global_mean.shape[0]), dtype = np.float32)
        for i, r in enumerate(idx):
            key = tuple(META[k].iloc[r] for k in self.keys)
            key = key[0] if len(key) == 1 else key
            out[i] = self.table.get(key, self.global_mean)
        return out


def permute_embedding(E: np.ndarray, META: pd.DataFrame, rng,
                      max_tries: int = 200) -> np.ndarray:
    # shuffles which compound gets which embeddings but preserves the distro
    perts = META["perturbation"].values
    uniq = pd.unique(perts)
    for _ in range(max_tries):
        shuffled = rng.permutation(uniq)
        if not np.any(shuffled == uniq):
            break
    else:
        raise RuntimeError("couldn't  find a derangement")

    fake_mapping = dict(zip(uniq, shuffled))
    first_row = {p: np.where(perts == p)[0][0] for p in uniq}
    shuffled_embeddings = np.vstack([E[first_row[fake_mapping[p]]] for p in perts])
    return shuffled_embeddings


def add_dose_feature(E: np.ndarray, 
                     dose: np.ndarray, ref: Optional[np.ndarray] = None,
                     interaction: bool = True):
    # calculats d = log_10(d + 1)
    ld = np.log10(dose.astype(np.float64) + 1.0)[:, None]
    src = np.log10(np.asarray(ref, dtype = np.float64) + 1.0) if ref is not None else ld
    ld = ((ld - src.mean()) / (src.std() + 1e-8)).astype(E.dtype)
    parts = [E, ld] + ([E * ld] if interaction else [])
    embedding_plus_dose = np.concatenate(parts, axis = 1)
    return embedding_plus_dose


# =========
# SECTION 5
# metrics helper functions
# 1. r2 pooled - overall pred accuracy across responses and genes
# 2. r2 vs baseline - improvement of model over cell line + dosage baseline 
# 3. row correlation - pearson correlation for each individual perturbation response 
# 4. test cipher hypotehssi to check if the most predictable response can be predicted w 
#       the first few high variance modes 
# =========+=
def r2_pooled(pred, true) -> float:
    ss_res = ((pred - true) ** 2).sum()
    ss_tot = ((true - true.mean(0, keepdims = True)) ** 2).sum()
    return float(1.0 - ss_res / ss_tot)


def r2_vs_baseline(pred, true, base) -> float:
    ss_res = ((pred - true) ** 2).sum()
    ss_base = ((base - true) ** 2).sum()
    contribution = float(1.0 - ss_res / ss_base)
    return contribution 


def row_correlation(a, b) -> np.ndarray:
    a = a - a.mean(1, keepdims = True)
    b = b - b.mean(1, keepdims = True)
    num = (a * b).sum(1)
    den = np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1)) + 1e-12
    return num / den


# eval
def evaluate(pred, true, base) -> Dict[str, float]:
    rc = row_correlation(pred, true)
    rr = row_correlation(pred - base, true - base)
    return dict(
        r2_pooled = r2_pooled(pred, true),
        r2_vs_baseline = r2_vs_baseline(pred, true, base),
        pearson_row_mean = float(rc.mean()),
        pearson_row_median = float(np.median(rc)),
        pearson_resid_mean =  float(rr.mean()),
        pearson_resid_median = float(np.median(rr)),
    )


def mode_resolved_r2(pred, true, basis: FluctuationBasis, ks = (1, 3, 10, 30, 100, 200, 500)):
    # check if performance here verifies the original paper's claim too or nah 
    P, T = basis.to_modes(pred), basis.to_modes(true)
    out, lo = {}, 0
    for k in ks:
        hi = min(k, P.shape[1])
        if hi <= lo:
            continue
        num = ((P[:, lo:hi] - T[:, lo:hi]) ** 2).sum()
        den = ((T[:, lo:hi] - T[:, lo:hi].mean(0, keepdims=True)) ** 2).sum() + 1e-12
        out[f"r2_modes_{lo}_{hi}"] = float(1 - num / den)
        lo = hi
    return out

# =========
# SECTION 6
# claude: nested CV 
# claude: 
    # Outer loop: unseen-compound folds, used only for reporting.
    # Inner loop: unseen-compound folds *within* the outer training set, used to
    # pick (shrinkage, lambda, gamma). Nothing about the outer test set touches
    # model selection -- which is the bug that made the old NN numbers
    # uninterpretable, since it early-stopped on the test fold and therefore
    # reported the best-of-150 test scores for real and shuffled embeddings alike.
# =========
def _fit_predict_cellwise(E, DX, META, train_mask, test_mask, bases,
                          shrinkage, lambdas, gammas, cfg: Config,
                          base_tr = None, base_te = None, in_modes: bool = False):
    cl_all = META["cell_line"].values
    n_te = int(test_mask.sum())
    preds = {(l, g): np.zeros((n_te, DX.shape[1]), np.float32)
             for g in gammas for l in lambdas}
    te_pos = np.zeros(len(META), int)
    te_pos[test_mask] = np.arange(n_te)

    target = DX if base_tr is None else DX - _scatter(base_tr, base_te, train_mask, test_mask, DX.shape)

    for cl in np.unique(cl_all[test_mask]):
        tr = train_mask & (cl_all == cl)
        te = test_mask & (cl_all == cl)
        if tr.sum() < 5 or te.sum() == 0:
            continue
        mr = ModeRidge(bases[cl], shrinkage).fit(E[tr], target[tr])
        rows = te_pos[te]
        for lam, gamma, p in mr.predict_grid(E[te], lambdas, gammas, in_modes=in_modes):
            preds[(lam, gamma)][rows] = p
    return preds


def _scatter(base_tr, base_te, train_mask, test_mask, shape):
    out = np.zeros(shape, np.float32)
    out[train_mask] = base_tr
    out[test_mask] = base_te
    return out


def run_cv(E_named: Dict[str, np.ndarray], DX: np.ndarray, META: pd.DataFrame,
           bases: Dict[str, FluctuationBasis], cfg: Config,
           select_on: str = "r2_vs_baseline",
           perm_controls: Sequence[str] = (),
           n_perm: int = 1) -> pd.DataFrame:
    from sklearn.model_selection import GroupKFold

    assert isinstance(META.index, pd.RangeIndex), "META must be reset_index'd; rows align with DX"
    rng = np.random.default_rng(cfg.random_state)
    lambdas, gammas = list(cfg.lambda_grid), list(cfg.gamma_grid)
    records = []

    for fold in sorted(META["fold"].unique()):
        test_mask = (META["fold"] == fold).values
        train_mask = ~test_mask
        tr_idx = np.where(train_mask)[0]

        sm = StratumMean().fit(DX, META, train_mask)
        base_te = sm.predict(META, np.where(test_mask)[0])
        base_tr = sm.predict(META, tr_idx)
        true_te = DX[test_mask]

        records.append(dict(fold=fold, condition="BASELINE__stratum_mean",
                            shrinkage=np.nan, lam=np.nan, gamma=np.nan,
                            **evaluate(base_te, true_te, base_te)))

        # inner folds, grouped by compound within the outer training set
        inner = GroupKFold(n_splits=cfg.n_inner_folds)
        inner_splits = list(inner.split(tr_idx, groups=META["perturbation"].values[tr_idx]))

        conditions = dict(E_named)
        for k in perm_controls:
            if k not in E_named:
                continue
            for j in range(n_perm):
                conditions[f"PERM{j:02d}__{k}"] = permute_embedding(E_named[k], META, rng)

        for name, E_full in conditions.items():
            E = add_dose_feature(E_full, META["dose_value"].values,
                                 ref=META["dose_value"].values[train_mask],
                                 interaction=cfg.dose_interaction)

            # ---- inner selection ----
            scores = {}
            for sh in cfg.shrinkage_grid:
                agg = {(l, g): [] for g in gammas for l in lambdas}
                for itr_rel, ite_rel in inner_splits:
                    itr = np.zeros(len(META), bool); itr[tr_idx[itr_rel]] = True
                    ite = np.zeros(len(META), bool); ite[tr_idx[ite_rel]] = True
                    sm_i = StratumMean().fit(DX, META, itr)
                    b_tr = sm_i.predict(META, np.where(itr)[0])
                    b_te = sm_i.predict(META, np.where(ite)[0])
                    fast = select_on in ("r2_vs_baseline", "r2_pooled")
                    preds = _fit_predict_cellwise(
                        E, DX, META, itr, ite, bases, sh, lambdas, gammas, cfg,
                        base_tr=b_tr if cfg.residual_formulation else None,
                        base_te=b_te if cfg.residual_formulation else None,
                        in_modes=fast)
                    if fast:
                        cl_i = META["cell_line"].values[ite]
                        tgt = np.zeros_like(preds[next(iter(preds))])
                        resid = DX[ite] - b_te if cfg.residual_formulation else DX[ite]
                        for cl in np.unique(cl_i):
                            m = cl_i == cl
                            tgt[m] = bases[cl].to_modes(resid[m])
                        ss_base = (resid ** 2).sum()
                        for key, p in preds.items():
                            agg[key].append(float(1 - ((p - tgt) ** 2).sum() / ss_base))
                    else:
                        for key, p in preds.items():
                            p_full = p + b_te if cfg.residual_formulation else p
                            agg[key].append(evaluate(p_full, DX[ite], b_te)[select_on])
                for key, vals in agg.items():
                    scores[(sh,) + key] = float(np.mean(vals))
            best_sh, best_lam, best_gamma = max(scores, key=scores.get)

            preds = _fit_predict_cellwise(
                E, DX, META, train_mask, test_mask, bases, best_sh,
                [best_lam], [best_gamma], cfg,
                base_tr=base_tr if cfg.residual_formulation else None,
                base_te=base_te if cfg.residual_formulation else None)
            pred = preds[(best_lam, best_gamma)]
            if cfg.residual_formulation:
                pred = pred + base_te

            rec = dict(fold=fold, condition=name, shrinkage=best_sh,
                       lam=best_lam, gamma=best_gamma,
                       inner_score=scores[(best_sh, best_lam, best_gamma)],
                       **evaluate(pred, true_te, base_te))
            records.append(rec)
            print(f"fold {fold} {name:28s} gamma={best_gamma:.2f} "
                  f"lam={best_lam:.3g} sh={best_sh:.2f} "
                  f"dR2={rec['r2_vs_baseline']:+.4f} "
                  f"r_resid={rec['pearson_resid_mean']:+.3f}")

    return pd.DataFrame(records)


def leaderboard(res: pd.DataFrame) -> pd.DataFrame:
    cols = ["r2_pooled", "r2_vs_baseline", "pearson_row_mean", "pearson_resid_mean"]
    lb = res.groupby("condition")[cols].agg(["mean", "std"])
    lb.columns = [f"{a}_{b}" for a, b in lb.columns]
    return lb.sort_values("r2_vs_baseline_mean", ascending=False).reset_index()


# # MAIN!! 