"""Score saved CellFlow predictions separately from training/prediction.

Predictions are read from backbone_eval/config.RESULTS/predictions/{spec}/fold_{fold}.
The observed held-out data and fold-specific PCA basis are reconstructed from the
cached preprocessed sciPlex3 dataset, so scoring does not require retraining CellFlow.

Run: 

PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/cellflow/bin/python -u backbone_eval/score.py     --specs [NAME]     --folds 0 1 2 3 4

CUDA_VISIBLE_DEVICES=1 PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/cellflow/bin/python -u backbone_eval/score.py     --specs sand/z   --folds 0 1 2 3 4

# to see all the scores calculated so far 
PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/core/bin/python -c '
import pandas as pd
pd.set_option("display.width", 240)
df = pd.read_parquet("data/results/backbone/backbone_metrics.parquet")
print(len(df), "rows\n")
print(df.to_string(index=False))
print("\nper spec:")
print(df.groupby("spec")["test_mmd"].agg(["mean","std","count"]).round(4).to_string())'

test_r_squared


# for each condition, use the fold's control cells as the "prediction"
PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/core/bin/python -c '
pair = ["fusion/pca_concat50/sand_x_combined_semantic", "ecfp_raw/counts"]
w = (cond[cond.spec.isin(pair)]
     .pivot_table(index=["fold", "condition_id"], columns="spec",
                  values="test_r_squared").dropna())
d = w[pair[0]] - w[pair[1]]
print(f"n={len(d)}  mean {d.mean():+.5f}  ci95 ±{1.96*d.std()/len(d)**0.5:.5f}")
print(f"wins {(d>0).mean():.1%}  wilcoxon p={stats.wilcoxon(w[pair[0]], w[pair[1]])[1]:.4f}")'


ECFP RAW E DISTANCE = 9.6944798
10.301250
15.916899
5.456076
9.783479
7.014695

PCA CONCAT E DISTANCE = 9.077877 = lower e distance!!! yayy 
8.956654
14.794541
6.378985
7.930052
7.329153
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from backbone_eval import config as cfg
from backbone_eval.data import compound_folds, load_adata, split

RUNS_PATH = cfg.RESULTS / "backbone_runs.parquet"
METRICS_PATH = cfg.RESULTS / "backbone_metrics.parquet"
CONDITION_METRICS_PATH = cfg.RESULTS / "backbone_condition_metrics.parquet"
PREDICTIONS_ROOT = cfg.RESULTS / "predictions"


def score(
    specs: list[str],
    folds: tuple[int, ...],
    force: bool = False,
) -> pd.DataFrame:
    from cellflow.metrics import compute_mean_metrics, compute_metrics

    print("[score] loading cached data", flush=True)
    adata = load_adata()
    fold_map = compound_folds(adata)

    all_aggregate_rows: list[dict] = []
    all_condition_rows: list[dict] = []

    for fold in folds:
        print(f"[score] FOLD {fold}", flush=True)
        t0 = time.perf_counter()
        _, test = split(adata, fold_map[fold])
        print(
            f"[score] split + PCA fold {fold}: "
            f"{time.perf_counter() - t0:.2f}s",
            flush=True,
        )

        perturbed = test[~test.obs["is_control"]].copy()
        pred_base = PREDICTIONS_ROOT

        for spec in specs:
            pred_dir = pred_base / spec / f"fold_{fold}"
            conditions_path = pred_dir / "conditions.parquet"
            metadata_path = pred_dir / "metadata.json"

            if not conditions_path.exists():
                raise FileNotFoundError(
                    f"No saved conditions for {spec} fold {fold}: {conditions_path}"
                )
            if not metadata_path.exists():
                raise FileNotFoundError(
                    f"No saved metadata for {spec} fold {fold}: {metadata_path}"
                )

            metadata = json.loads(metadata_path.read_text())
            conditions = pd.read_parquet(conditions_path)

            print(
                f"[score] {spec} fold {fold}: "
                f"{len(conditions)} conditions",
                flush=True,
            )

            per_condition: dict[str, dict] = {}
            condition_rows: list[dict] = []

            t_score = time.perf_counter()
            for j, row in conditions.iterrows():
                condition_id = row["condition_id"]
                pred_path = pred_dir / row["prediction_file"]

                if not pred_path.exists():
                    raise FileNotFoundError(
                        f"Missing prediction for {condition_id}: {pred_path}"
                    )

                predicted = np.load(pred_path, allow_pickle=False)

                mask = (
                    (perturbed.obs["perturbation"] == row["perturbation"])
                    & (perturbed.obs["log_dose"] == row["log_dose"])
                    & (perturbed.obs["cell_line"] == row["cell_line"])
                )
                observed = perturbed[mask.values].obsm["X_pca"]

                if len(observed) == 0:
                    print(
                        f"[score] {j + 1}/{len(conditions)} "
                        f"{condition_id}: no observed cells, skip",
                        flush=True,
                    )
                    continue

                t0 = time.perf_counter()
                raw_metrics = compute_metrics(observed, predicted)
                dt = time.perf_counter() - t0

                metrics = {k: float(v) for k, v in raw_metrics.items()}
                per_condition[condition_id] = metrics

                condition_rows.append(
                    {
                        "spec": spec,
                        "fold": fold,
                        "condition_id": condition_id,
                        "n_iterations": metadata["n_iterations"],
                        "predict_rtol": metadata["predict_rtol"],
                        "predict_atol": metadata["predict_atol"],
                        "n_observed": len(observed),
                        "n_predicted": len(predicted),
                        "metric_seconds": dt,
                        **{
                            f"test_{k}" if not k.startswith("test_") else k: v
                            for k, v in metrics.items()
                        },
                    }
                )

                print(
                    f"[score] {j + 1}/{len(conditions)} "
                    f"{condition_id} | {dt:.2f}s",
                    flush=True,
                )

            all_condition_rows.extend(condition_rows)

            aggregate = compute_mean_metrics(per_condition, prefix="test_")
            aggregate = {k: float(v) for k, v in aggregate.items()}

            aggregate_row = {
                "spec": spec,
                "fold": fold,
                "n_iterations": metadata["n_iterations"],
                "predict_rtol": metadata["predict_rtol"],
                "predict_atol": metadata["predict_atol"],
                "n_conditions": len(per_condition),
                "metric_seconds": round(time.perf_counter() - t_score, 1),
                **aggregate,
            }
            all_aggregate_rows.append(aggregate_row)

            print(
                f"[score] DONE {spec} fold {fold} | "
                f"{aggregate} | "
                f"metrics {aggregate_row['metric_seconds']:.1f}s",
                flush=True,
            )
            
            _upsert_rows(
                METRICS_PATH,
                pd.DataFrame([aggregate_row]),
                keys=["spec", "fold", "n_iterations", "predict_rtol", "predict_atol"],
                force=force,
            )

            _upsert_rows(
                CONDITION_METRICS_PATH,
                pd.DataFrame(condition_rows),
                keys=[
                    "spec",
                    "fold",
                    "condition_id",
                    "n_iterations",
                    "predict_rtol",
                    "predict_atol",
                ],
                force=force,
            )
    return pd.DataFrame(all_aggregate_rows)


def _upsert_rows(
    path: Path,
    new: pd.DataFrame,
    keys: list[str],
    force: bool,
) -> None:
    if new.empty:
        return

    if path.exists() and not force:
        old = pd.read_parquet(path)
    elif path.exists():
        old = pd.read_parquet(path)
    else:
        old = pd.DataFrame()

    if not old.empty:
        # Remove rows matching the keys supplied by the new rows.
        for _, row in new[keys].drop_duplicates().iterrows():
            mask = np.ones(len(old), dtype=bool)
            for key in keys:
                mask &= old[key].eq(row[key]).to_numpy()
            old = old.loc[~mask]

    out = pd.concat([old, new], ignore_index=True)
    out.to_parquet(path, index=False)


def summary() -> pd.DataFrame:
    if not METRICS_PATH.exists():
        return pd.DataFrame()
    df = pd.read_parquet(METRICS_PATH)
    if df.empty or "test_r_squared" not in df.columns:
        return df

    out = (
        df.groupby("spec")["test_r_squared"]
        .agg(mean="mean", sd="std", n="count")
        .reset_index()
    )
    out["ci95"] = 1.96 * out["sd"] / np.sqrt(out["n"])

    if cfg.BASELINE in set(out["spec"]):
        base = float(out.loc[out["spec"] == cfg.BASELINE, "mean"].iloc[0])
        out["vs_baseline"] = out["mean"] - base

    return out.sort_values("mean", ascending=False)


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--specs", nargs="+", required=True)
    p.add_argument("--folds", nargs="+", type=int, required=True)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    df = score(args.specs, tuple(args.folds), force=args.force)
    print("\nPer-fold metrics:")
    print(df.to_string(index=False))
    print("\nSummary:")
    print(summary().to_string(index=False))