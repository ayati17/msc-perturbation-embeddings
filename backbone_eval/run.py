"""
Train CellFlow and save held-out predictions, without computing metrics bc it was taking too long together.

One training run is performed per (embedding, fold). 
Predictions are saved per condition
Scoring is handled separately by score.py.
"""

from __future__ import annotations
from contextlib import contextmanager
from functools import partial
import json
import shutil
import time
from pathlib import Path
import diffrax
import numpy as np
import pandas as pd

from backbone_eval import config as cfg
from backbone_eval.data import compound_folds, load_adata, split
from backbone_eval.reps import drug_reps


PREDICT_RTOL = 1e-2
PREDICT_ATOL = 1e-2
# set to None for the full fold.
MAX_PREDICTION_CONDITIONS: int | None = None
# nt | None = None

RUNS_PATH = cfg.RESULTS / "backbone_runs.parquet"
PREDICTIONS_ROOT = cfg.RESULTS / "predictions"


def run(
    specs: list[str],
    folds: tuple | None = None,
    n_iterations: int = cfg.N_ITERATIONS,
    force: bool = False,
    out_name: str = "backbone",
) -> pd.DataFrame:
    """
    Train and predict once per (spec, fold) for cellflow!
    """
    del out_name  

    done = _completed(RUNS_PATH) if (RUNS_PATH.exists() and not force) else set()

    print("START", flush=True)
    adata = load_adata()
    print("[run] finished load_adata()", flush=True)
    print("[run] starting compound_folds()", flush=True)
    fold_map = compound_folds(adata)
    print("[run] finished compound_folds()", flush=True)
    folds = folds if folds is not None else tuple(fold_map)

    for fold in folds:
        print(f"FOLD {fold}", flush=True)

        with timer(f"split + PCA fold {fold}"):
            train, test = split(adata, fold_map[fold])

        print(f"STARTING SPECS FOR FOLD {fold}", flush=True)
        names = sorted(
            set(adata.obs.loc[~adata.obs["is_control"], "perturbation"])
        )

        for spec in specs:
            run_key = (
                spec,
                fold,
                n_iterations,
                PREDICT_RTOL,
                PREDICT_ATOL,
            )

            if run_key in done:
                print(f"[skip] {spec} fold {fold}", flush=True)
                continue

            started = time.perf_counter()

            with timer(f"load reps {spec}"):
                reps = drug_reps(spec, names)

            pred_dir = PREDICTIONS_ROOT / spec / f"fold_{fold}"
            if force and pred_dir.exists():
                shutil.rmtree(pred_dir)
            pred_dir.mkdir(parents=True, exist_ok=True)

            metadata_path = pred_dir / "metadata.json"
            conditions_path = pred_dir / "conditions.parquet"

            _train_and_predict(
                train=train,
                test=test,
                reps=reps,
                n_iterations=n_iterations,
                pred_dir=pred_dir,
                metadata_path=metadata_path,
                conditions_path=conditions_path,
                spec=spec,
                fold=fold,
            )

            seconds = time.perf_counter() - started
            row = {
                "spec": spec,
                "fold": fold,
                "n_dims": len(next(iter(reps.values()))),
                "n_iterations": n_iterations,
                "predict_rtol": PREDICT_RTOL,
                "predict_atol": PREDICT_ATOL,
                "seconds": round(seconds, 1),
                "prediction_dir": str(pred_dir),
                "n_prediction_conditions": _count_prediction_files(pred_dir),
            }
            _upsert_run_row(RUNS_PATH, row)

            print(
                f"[done] {spec:34s} fold {fold} "
                f"({seconds:.0f}s) predictions saved to {pred_dir}",
                flush=True,
            )

    return load_results()


def _train_and_predict(
    train,
    test,
    reps: dict,
    n_iterations: int,
    pred_dir: Path,
    metadata_path: Path,
    conditions_path: Path,
    spec: str,
    fold: int,
) -> None:
    from cellflow.model import CellFlow
    from cellflow.utils import match_linear

    train.uns["drug_reps"] = reps
    test.uns["drug_reps"] = reps

    with timer("CellFlow init"):
        model = CellFlow(train, solver="otfm")

    with timer("prepare_data"):
        model.prepare_data(
            sample_rep="X_pca",
            control_key="is_control",
            perturbation_covariates={"drug": ("perturbation",)},
            perturbation_covariate_reps={"drug": "drug_reps"},
            sample_covariates=["log_dose"],
            split_covariates=["cell_line"],
            max_combination_length=1,
            null_value=0.0,
        )

    with timer("prepare_model"):
        model.prepare_model(
            match_fn=partial(match_linear, tau_a=1.0, tau_b=1.0),
            **cfg.MODEL_KWARGS,
        )

    with timer("TRAIN"):
        model.train(
            num_iterations=n_iterations,
            batch_size=cfg.BATCH_SIZE,
        )

    controls = test[test.obs["is_control"]].copy()
    perturbed = test[~test.obs["is_control"]].copy()

    conditions = _condition_table(perturbed)
    if MAX_PREDICTION_CONDITIONS is not None:
        conditions = conditions.head(MAX_PREDICTION_CONDITIONS).copy()

    conditions["prediction_file"] = [
        f"pred_{i:04d}.npy" for i in range(len(conditions))
    ]
    conditions.to_parquet(conditions_path, index=False)

    print(
        f"[predict] {len(conditions)} conditions | "
        f"{len(controls):,} control cells | "
        f"rtol={PREDICT_RTOL:g} atol={PREDICT_ATOL:g}",
        flush=True,
    )

    metadata = {
        "spec": spec,
        "fold": fold,
        "n_iterations": n_iterations,
        "predict_rtol": PREDICT_RTOL,
        "predict_atol": PREDICT_ATOL,
        "n_conditions": len(conditions),
        "n_control_cells": len(controls),
        "n_dims": len(next(iter(reps.values()))),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2))

    with timer("PREDICT"):
        for i, condition_row in conditions.iterrows():
            condition = conditions.loc[[i]]
            condition_id = condition_row["condition_id"]
            pred_path = pred_dir / condition_row["prediction_file"]

            if pred_path.exists():
                print(
                    f"[predict] {i + 1}/{len(conditions)} {condition_id} [saved; skip]",
                    flush=True,
                )
                continue

            print(
                f"[predict] {i + 1}/{len(conditions)} {condition_id}",
                flush=True,
            )

            t0 = time.perf_counter()
            pred = model.predict(
                controls,
                covariate_data=condition,
                sample_rep="X_pca",
                condition_id_key="condition_id",
                stepsize_controller=diffrax.PIDController(
                    rtol=PREDICT_RTOL,
                    atol=PREDICT_ATOL,
                ),
            )
            dt = time.perf_counter() - t0

            # model.predict() returns a dict keyed by condition_id.
            arr = np.asarray(pred[condition_id])
            np.save(pred_path, arr)

            print(
                f"[predict] {condition_id} finished in {dt:.2f}s "
                f"-> {pred_path.name}",
                flush=True,
            )


def _condition_table(perturbed) -> pd.DataFrame:
    """ One row per (compound, dose, cell line) to predict."""
    conditions = (
        perturbed.obs[["perturbation", "log_dose", "cell_line"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    conditions["is_control"] = False
    conditions["condition_id"] = (
        conditions["perturbation"].astype(str)
        + "__"
        + conditions["log_dose"].astype(str)
        + "__"
        + conditions["cell_line"].astype(str)
    )
    return conditions


def _upsert_run_row(path: Path, row: dict) -> None:
    df_new = pd.DataFrame([row])
    if path.exists():
        df = pd.read_parquet(path)
        keys = ["spec", "fold", "n_iterations", "predict_rtol", "predict_atol"]
        mask = np.ones(len(df), dtype=bool)
        for key in keys:
            mask &= df[key].eq(row[key]).to_numpy()
        df = df.loc[~mask]
        df_new = pd.concat([df, df_new], ignore_index=True)
    df_new.to_parquet(path, index=False)


def _completed(path: Path) -> set:
    if not path.exists():
        return set()
    df = pd.read_parquet(path)
    keys = ["spec", "fold", "n_iterations", "predict_rtol", "predict_atol"]
    return set(map(tuple, df[keys].drop_duplicates().values))


def _count_prediction_files(pred_dir: Path) -> int:
    return len(list(pred_dir.glob("pred_*.npy")))


def load_results() -> pd.DataFrame:
    return pd.read_parquet(RUNS_PATH) if RUNS_PATH.exists() else pd.DataFrame()


def summary(df: pd.DataFrame | None = None):
    """Summarize completed train+predict runs."""
    df = load_results() if df is None else df
    if df.empty:
        return df
    return df.sort_values(["spec", "fold", "n_iterations"])


@contextmanager
def timer(name):
    t0 = time.perf_counter()
    yield
    dt = time.perf_counter() - t0
    print(f"[TIME] {name:25s} {dt:10.2f}s", flush=True)