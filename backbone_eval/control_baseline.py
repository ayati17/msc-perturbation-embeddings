"""
Dummy prediction round

    backbone_eval/control_baseline.py
    CUDA_VISIBLE_DEVICES=2 PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/cellflow/bin/python -u \
        backbone_eval/control_baseline.py --folds 0 1 2 3 4
"""

from __future__ import annotations

import argparse
import time
import numpy as np
import pandas as pd

from backbone_eval import config as cfg
from backbone_eval.data import compound_folds, load_adata, split

OUT = cfg.RESULTS / "backbone_control_baseline.parquet"


def run(folds: tuple[int, ...], seed: int = 17) -> pd.DataFrame:
    from cellflow.metrics import compute_metrics

    adata = load_adata()
    fold_map = compound_folds(adata)
    rng = np.random.default_rng(seed)

    rows = []
    for fold in folds:
        started = time.perf_counter()
        _, test = split(adata, fold_map[fold])

        controls = test[test.obs["is_control"]].obsm["X_pca"]
        perturbed = test[~test.obs["is_control"]].copy()

        conditions = (perturbed.obs[["perturbation", "log_dose", "cell_line"]]
                      .drop_duplicates().reset_index(drop=True))
        print(f"[baseline] fold {fold}: {len(conditions)} conditions, "
              f"{len(controls)} control cells", flush=True)

        for _, condition in conditions.iterrows():
            mask = ((perturbed.obs["perturbation"] == condition.perturbation)
                    & (perturbed.obs["log_dose"] == condition.log_dose)
                    & (perturbed.obs["cell_line"] == condition.cell_line))
            observed = perturbed[mask.values].obsm["X_pca"]
            if len(observed) == 0:
                continue

            # the do-nothing prediction: control cells of the same cell line,
            # sampled to match the observed count
            same_line = (test.obs["is_control"]
                         & (test.obs["cell_line"] == condition.cell_line)).to_numpy()
            pool = test[same_line].obsm["X_pca"]
            picked = pool[rng.choice(len(pool), size=len(observed),
                                     replace=len(pool) < len(observed))]

            metrics = {f"test_{k}": float(v)
                       for k, v in compute_metrics(observed, picked).items()}
            rows.append(dict(
                spec="CONTROL__do_nothing", fold=fold,
                condition_id=f"{condition.perturbation}__{condition.log_dose}"
                             f"__{condition.cell_line}",
                perturbation=condition.perturbation, cell_line=condition.cell_line,
                log_dose=condition.log_dose, n_observed=len(observed), **metrics))

        print(f"[baseline] fold {fold} done in "
              f"{time.perf_counter() - started:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_parquet(OUT, index=False)
    print(f"\nsaved {len(df)} rows -> {OUT}")
    return df


def report(df: pd.DataFrame | None = None) -> None:
    df = pd.read_parquet(OUT) if df is None else df
    metrics = [c for c in df.columns if c.startswith("test_")]

    print("\ndo-nothing baseline, per fold:")
    print(df.groupby("fold")[metrics].mean().round(4).to_string())

    print("\noverall:")
    print(df[metrics].mean().round(4).to_string())

    metrics_path = cfg.RESULTS / "backbone_metrics.parquet"
    if metrics_path.exists() and "test_r_squared" in df.columns:
        trained = pd.read_parquet(metrics_path)
        baseline = df.groupby("fold")["test_r_squared"].mean().rename("do_nothing")
        merged = (trained.groupby(["spec", "fold"])["test_r_squared"].mean()
                         .reset_index().merge(baseline, on="fold"))
        merged["gain"] = merged.test_r_squared - merged.do_nothing
        print("\ntrained models vs do-nothing, paired by fold:")
        print(merged.groupby("spec")[["test_r_squared", "do_nothing", "gain"]]
                    .mean().round(4).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    args = parser.parse_args()
    report(run(tuple(args.folds)))