from __future__ import annotations

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .bridge import load_from_artifacts
from .cipher_lr import (
    Config, ModeRidge, StratumMean, add_dose_feature, assign_outer_folds,
    build_bases, build_response_matrix, evaluate, leaderboard,
    mode_resolved_r2, observed_force, run_cv,
)

DATASET = "sciplex3"

def spec_key(spec: str) -> str:
    return spec.replace("/", "__")

def cache_file(cfg: Config) -> Path:
    return Path(cfg.work_dir) / "cipher_cache.pkl"


def make_config() -> Config:
    return Config()

def known_compound_names(dataset: str = DATASET) -> set[str]:

    from chemembed.config import COMPOUNDS

    table = pd.read_csv(COMPOUNDS / f"{dataset}.csv")
    return set(table["dataset_drug_name"].astype(str).str.strip())


def load_everything(cfg: Config, which, use_cache: bool = True):
    # return DX, META, bases, E_named, genes
    which = list(which) if which else list(cfg.embeddings)
    if not which:
        sys.exit("no embeddings requested -- pass --use SPEC1,SPEC2")

    cache = cache_file(cfg)
    if use_cache and cache.exists():
        print(f"loading cache {cache}")
        blob = pickle.load(open(cache, "rb"))
    else:
        import scanpy as sc

        t0 = time.time()
        adata = sc.read_h5ad(cfg.adata_path)
        print(f"h5ad loaded in {time.time() - t0:.0f}s: {adata.shape}")

        DX, META, X, obs, genes, is_control = build_response_matrix(adata, cfg)

        known = known_compound_names(cfg.dataset)
        keep = np.array([p in known for p in META["perturbation"].values])
        if not keep.all():
            lost = sorted(set(META["perturbation"].values[~keep]))
            print(f"\ndropping {int((~keep).sum())} rows, {len(lost)} compounds absent "
                  f"from {cfg.dataset}.csv: {lost[:8]}{' ...' if len(lost) > 8 else ''}")
        META = META.loc[keep].reset_index(drop=True)
        DX = DX[keep]

        META = assign_outer_folds(META, cfg)
        bases = build_bases(X, obs, is_control, sorted(META["cell_line"].unique()))

        blob = dict(DX=DX, META=META, bases=bases, genes=genes)
        with open(cache, "wb") as fh:
            pickle.dump(blob, fh)
        print(f"cached to {cache}")

    META = blob["META"]
    E_named = load_from_artifacts(which, META["perturbation"].values,
                                  dataset=cfg.dataset)
    missing = [s for s in which if spec_key(s) not in E_named]
    if missing:
        sys.exit(f"no artifact for: {missing}\nbuild them first, or drop them from --use")

    return blob["DX"], META, blob["bases"], E_named, blob["genes"]


def stage_diagnose(cfg: Config, which):
    DX, META, bases, E_named, genes = load_everything(cfg, which)

    print("\n=== response matrix ===")
    print(f"DX {DX.shape}   compounds={META['perturbation'].nunique()}   "
          f"doses={sorted(META['dose_value'].unique())}")
    print(META.groupby("cell_line")["n_cells"].describe()[["count", "min", "50%", "max"]])
    thin = (META["n_cells"] < 20).mean()
    print(f"fraction of rows built from <20 cells: {thin:.1%}")
    if thin > 0.25:
        print("  -> consider raising cfg.min_cells_per_group; these rows carry "
              "several times the sampling noise of the rest but equal weight")

    print("\n=== fold balance (unseen-compound split) ===")
    print(META.groupby("fold").agg(rows=("perturbation", "size"),
                                   compounds=("perturbation", "nunique")))

    print("\n=== fluctuation structure ===")
    for cl, b in bases.items():
        print(f"{cl}: n_ctrl={b.n_cells}  G={b.G}  rank<={min(b.n_cells - 1, b.G)}")
        for d in cfg.shrinkage_grid:
            s = b.shrink_eigenvals(d)
            print(f"    delta={d:<5} eff_rank={b.effective_rank(d):7.1f}  "
                  f"cond={s.max() / s.min():.3g}")
        print(f"    cumulative variance: {b.variance_dict((1, 3, 10, 30, 100))}")
        ratio = b.n_cells / b.G
        if b.n_cells <= b.G:
            print(f"    WARNING: {b.n_cells} control cells for {b.G} genes -- Sigma is "
                  f"rank-deficient (rank <= {b.n_cells - 1}); shrinkage is the only "
                  f"thing making Sigma^-gamma defined, so results describe the "
                  f"shrinkage target as much as the data")
        elif ratio < 10:
            print(f"    note: n/G = {ratio:.1f}. Sigma is full rank, but Marchenko-"
                  f"Pastur still spreads the sample eigenvalues well beyond the truth "
                  f"at this ratio -- shrinkage is about spectrum accuracy, not "
                  f"invertibility")

    print("\n=== embeddings ===")
    for k, E in E_named.items():
        n_uniq = len(np.unique(E, axis=0))
        flag = "" if n_uniq == META["perturbation"].nunique() else \
            f"  <- expected {META['perturbation'].nunique()}; salt pairs share a compound_id"
        print(f"{k:28s} shape={E.shape}  finite={np.isfinite(E).all()}  "
              f"std={E.std():.4g}  n_unique_rows={n_uniq}{flag}")

    print("\n=== response scale ===")
    print(f"|DX| mean abs = {np.abs(DX).mean():.4g}, "
          f"per-row norm median = {np.median(np.linalg.norm(DX, axis=1)):.4g}")
    sm = StratumMean().fit(DX, META, np.ones(len(META), bool))
    base = sm.predict(META, np.arange(len(META)))
    frac = 1 - ((DX - base) ** 2).sum() / ((DX - DX.mean(0)) ** 2).sum()
    print(f"variance explained by (cell line, dose) alone: {frac:.3f}")
    print("  -> this is the share the chemistry model does NOT have to earn; "
          "everything is scored against it")

    n_fits = (cfg.n_outer_folds * cfg.n_inner_folds * len(cfg.shrinkage_grid)
              * (len(E_named) + 1))
    print(f"\nestimated SVD fits for `run`: ~{n_fits} x n_cell_lines")

def stage_run(cfg: Config, which, n_perm: int):
    DX, META, bases, E_named, genes = load_everything(cfg, which)
    out = Path(cfg.work_dir)

    t0 = time.time()
    perm_on = next(iter(E_named))
    res = run_cv(E_named, DX, META, bases, cfg,
                 perm_controls=(perm_on,), n_perm=n_perm)
    print(f"\ntotal {time.time() - t0:.0f}s")

    res.to_csv(out / "cv_results.csv", index=False)
    lb = leaderboard(res)
    lb.to_csv(out / "leaderboard.csv", index=False)
    (out / "config.txt").write_text(repr(cfg))

    print("\n=== leaderboard ===")
    print(lb.to_string(index=False))

    print("\n=== selected hyperparameters ===")
    sel = res.dropna(subset=["gamma"])
    print(sel.groupby("condition")[["gamma", "lam", "shrinkage"]]
          .agg(["median", "min", "max"]).to_string())

    perm = res[res["condition"].str.startswith("PERM")]
    real = res[~res["condition"].str.startswith(("PERM", "BASELINE"))]
    if len(perm):
        print(f"\nnull distribution over {n_perm} permutation(s) x "
              f"{cfg.n_outer_folds} folds (shuffled on {perm_on}):")
        print(f"  dR2       null {perm['r2_vs_baseline'].mean():+.4f} "
              f"+/- {perm['r2_vs_baseline'].std():.4f}   "
              f"(95th pct {np.percentile(perm['r2_vs_baseline'], 95):+.4f}, "
              f"max {perm['r2_vs_baseline'].max():+.4f})")
        print(f"  r_resid   null {perm['pearson_resid_mean'].mean():+.4f} "
              f"+/- {perm['pearson_resid_mean'].std():.4f}")
        for cond, g in real.groupby("condition"):
            z = ((g["r2_vs_baseline"].mean() - perm["r2_vs_baseline"].mean())
                 / (perm["r2_vs_baseline"].std() + 1e-12))
            print(f"  {cond:28s} dR2={g['r2_vs_baseline'].mean():+.4f}  z={z:+.2f}")
        # print("  (the nulls share folds, so treat z as a heuristic -- what matters is "
        #       "whether the real score clears every null draw)")
    return res, lb


# ---------------------------------------------------------------------
# stage 3: inspect
# ---------------------------------------------------------------------

def stage_inspect(cfg: Config, which, condition: str, fold: int):
    DX, META, bases, E_named, genes = load_everything(cfg, which)
    if condition not in E_named:
        sys.exit(f"--condition {condition} not loaded; available: {sorted(E_named)}")

    res = pd.read_csv(Path(cfg.work_dir) / "cv_results.csv")
    hit = res[(res["condition"] == condition) & (res["fold"] == fold)]
    if hit.empty:
        sys.exit(f"no row for {condition} fold {fold} in cv_results.csv; "
                 f"run the `run` stage first")
    row = hit.iloc[0]
    sh, lam, gamma = float(row["shrinkage"]), float(row["lam"]), float(row["gamma"])
    print(f"{condition} fold {fold}: delta={sh} lambda={lam:.4g} gamma={gamma}")

    test = (META["fold"] == fold).values
    train = ~test
    cl_all = META["cell_line"].values
    sm = StratumMean().fit(DX, META, train)
    b_tr = sm.predict(META, np.where(train)[0])
    b_te = sm.predict(META, np.where(test)[0])

    E = add_dose_feature(E_named[condition], META["dose_value"].values,
                         ref=META["dose_value"].values[train],
                         interaction=cfg.dose_interaction)
    pred = np.zeros((test.sum(), DX.shape[1]), np.float32)
    te_pos = np.zeros(len(META), int); te_pos[test] = np.arange(test.sum())
    tr_pos = np.zeros(len(META), int); tr_pos[train] = np.arange(train.sum())
    forces = {}
    for cl in np.unique(cl_all[test]):
        tr, te = train & (cl_all == cl), test & (cl_all == cl)
        target = DX[tr] - b_tr[tr_pos[tr]] if cfg.residual_formulation else DX[tr]
        mr = ModeRidge(bases[cl], sh).fit(E[tr], target)
        pred[te_pos[te]] = mr.predict(E[te], lam, gamma)
        forces[cl] = mr.implied_force(E[te], lam, gamma)
    if cfg.residual_formulation:
        pred = pred + b_te

    true = DX[test]
    print("\nheld-out metrics:",
          {k: round(v, 4) for k, v in evaluate(pred, true, b_te).items()})

    print("\n=== where the model predicts, by mode block ===")
    cl0 = sorted(bases)[0]
    sel = cl_all[test] == cl0
    print(f"(cell line {cl0}, {sel.sum()} held-out rows)")
    for k, v in mode_resolved_r2(pred[sel], true[sel], bases[cl0],
                                 ks=(1, 3, 10, 30, 100, 300, 1000)).items():
        print(f"  {k:18s} {v:+.4f}")
    print("  -> CIPHER's claim is that responses live in the first few blocks. "
          "If r2 is flat across blocks, the response is not low-dimensional here.")

    print("\n=== implied force: is u sparse or dense? ===")
    for cl, U in forces.items():
        a = np.abs(U)
        frac = np.sort(a, axis=1)[:, ::-1][:, :50].sum(1) / (a.sum(1) + 1e-12)
        p = a / (a.sum(1, keepdims=True) + 1e-12)
        pr = np.exp(-(p * np.log(p + 1e-300)).sum(1))
        print(f"  {cl}: top-50 genes hold {frac.mean():.1%} of |u|; "
              f"participation ratio {pr.mean():.0f} / {DX.shape[1]} genes")

    U = observed_force(bases[cl0], DX[cl_all == cl0], sh)
    top = np.argsort(-np.abs(U).mean(0))[:20]
    print(f"\ngenes most consistently forced across all compounds ({cl0}):")
    print(", ".join(str(g) for g in np.asarray(genes)[top]))
    print("  -> if these are a generic stress/proliferation set rather than "
          "drug-specific targets, u is mostly shared and hard to predict from structure")


# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="CIPHER linear-response pipeline")
    ap.add_argument("stage", choices=["diagnose", "run", "inspect"])
    ap.add_argument("--use", required=True,
                    help="comma-separated chemembed specs, e.g. 'A/iupac,BASE/ecfp_raw'")
    ap.add_argument("--adata", default=None)
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--dataset", default=None, help=f"compound table name (default {DATASET})")
    ap.add_argument("--n-perm", type=int, default=5)
    ap.add_argument("--condition", default=None,
                    help="which loaded embedding to inspect (defaults to the first)")
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--no-residual", action="store_true")
    ap.add_argument("--no-interaction", action="store_true",
                    help="drop the embedding x log-dose interaction term")
    ap.add_argument("--fresh", action="store_true", help="rebuild the cache")
    args = ap.parse_args()

    cfg = make_config()
    if args.adata:
        cfg.adata_path = args.adata
    if args.work_dir:
        cfg.work_dir = args.work_dir
    if args.dataset:
        cfg.dataset = args.dataset
    if args.no_residual:
        cfg.residual_formulation = False
    if args.no_interaction:
        cfg.dose_interaction = False

    specs = [x.strip() for x in args.use.split(",") if x.strip()]
    cfg.embeddings = tuple(specs)

    Path(cfg.work_dir).mkdir(parents=True, exist_ok=True)
    if args.fresh and cache_file(cfg).exists():
        cache_file(cfg).unlink()

    if args.stage == "diagnose":
        stage_diagnose(cfg, specs)
    elif args.stage == "run":
        stage_run(cfg, specs, args.n_perm)
    else:
        condition = args.condition or spec_key(specs[0])
        stage_inspect(cfg, specs, condition, args.fold)


if __name__ == "__main__":
    main()