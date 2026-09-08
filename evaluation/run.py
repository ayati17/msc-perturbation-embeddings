"""
Pipeline: Basically builds (spec x regression model x cell line x split x fold) -> 1 row.
rows are appended to a parquet shard as each spec finishes
"""

from __future__ import annotations
import time
from itertools import product
import numpy as np
import pandas as pd
from evaluation import config as cfg
from evaluation.data import evaluate, load_xy, make_folds
from evaluation.models import build
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category = ConvergenceWarning)

def run(specs, models, cell_lines = cfg.CELL_LINES, dose = cfg.PRIMARY_DOSE,
        splits = cfg.SPLITS, force = False, out_name = "metrics",
        multi_dose: bool = False) -> pd.DataFrame:

    path = cfg.RESULTS /f"{out_name}.parquet"
    done = completed(path, splits, cfg.N_FOLDS) if (path.exists() and not force) else set()
    
    for spec, cell_line in product(specs, cell_lines):
        try:
            if multi_dose:
                from evaluation.data import load_xy_doses
                X, Y, ids, control_mean = load_xy_doses(spec, cell_line)
                dose_label = "all"
            else:
                X, Y, ids, control_mean = load_xy(spec, cell_line, dose)
                dose_label = dose
        except (FileNotFoundError, KeyError) as e:
            print(f"[skip] {spec} @ {cell_line}: {e}")
            continue

        print(f"\n {spec} @ {cell_line}  X{X.shape} Y{Y.shape} ")

        for model_name in models:
            if (spec, model_name, cell_line, dose_label) in done:
                print(f"{model_name:10s} already done")
                continue

            started = time.perf_counter()
            new = run_one(spec, model_name, cell_line, dose_label, X, Y, ids, splits, control_mean = control_mean)
            score = np.mean([r["value"] for r in new if r["metric"] == "pearson_top"])
            print(f"  {model_name:10s} pearson_top = {score:+.3f}  "
                  f"({time.perf_counter() - started:.1f}s)")
            append(path, new)       

    return load_results(out_name)

def _model_kwargs(spec: str, model_name: str) -> dict:
    """Per-spec model arguments. Block-structured models need to know where
    the blocks are, which only fusion artifacts can tell them."""
    kwargs = dict(cfg.SPEC_KWARGS.get(spec, {}))
    needs_block_dims = (
        model_name.startswith("gp_additive")
        or model_name.startswith("so_pls")
        or model_name.startswith("mb_pls")
    )
    if needs_block_dims and spec.startswith("fusion/"):
        import json
        from chemembed.artifacts import artifact_dir
        meta = json.loads((artifact_dir(spec) / "meta.json").read_text())
        kwargs["feature_dims"] = tuple(meta["block_dims"])
    
    if needs_block_dims and not spec.startswith("fusion/"):
        raise ValueError(
            f"{model_name} needs block boundaries but {spec} is not a fusion "
            f"artifact — only fusion/*/* records block_dims in meta.json")
    return kwargs

def run_one(spec, model_name, cell_line, dose, X, Y, ids, splits, control_mean = None) -> list[dict]:
    """
    One row per (split, fold, metric), with the model's hyperparameters and the
    config that produced them recorded as columns.
    """
    import dataclasses
    import datetime
    
    kwargs = _model_kwargs(spec, model_name)
    if "feature_dims" in kwargs:
        assert sum(kwargs["feature_dims"]) == X.shape[1], (
            f"{spec}/{model_name}: feature_dims sums to {sum(kwargs['feature_dims'])} "
            f"but X has {X.shape[1]} columns — block concat order/meta.json mismatch"
        )

    rows = []
    for split_kind in splits:
        for fold, (train, test) in enumerate(make_folds(ids, split_kind)):
            model = build(model_name, **kwargs)
            model.fit(X[train], Y[train])
            
            test_predicted = model.predict(X[test])
            assert test_predicted.shape == Y[test].shape, (
                f"{model_name} on {spec}: predicted {test_predicted.shape}, "
                f"expected {Y[test].shape}")
            
            extra = {}
            if hasattr(model, "block_importance_"):
                extra["block_importance"] = str(np.round(model.block_importance_, 3).tolist())

            metrics = evaluate(Y[test], test_predicted, control_mean)
            train_metrics = evaluate(Y[train], model.predict(X[train]), control_mean)
            metrics |= {f"train_{k}": v for k, v in train_metrics.items()}

            params = (dataclasses.asdict(model)
                      if dataclasses.is_dataclass(model) else {})
            params = {f"p_{k}": (str(v) if isinstance(v, tuple) else v)
                      for k, v in params.items()
                      if v is None or isinstance(v, (int, float, str, bool, tuple))}

            base = dict(
                spec=spec, model=model_name,
                cell_line=cell_line, dose=dose,
                split=split_kind, fold=fold,
                n_train=len(train), n_test=len(test), n_features=X.shape[1],
                cfg_n_hvg=cfg.N_HVG, cfg_time=cfg.TIME,
                cfg_top_n_degs=cfg.TOP_N_DEGS, cfg_n_folds=cfg.N_FOLDS,
                run_date=datetime.date.today().isoformat(),
                **params, **extra,
            )
            rows += [dict(base, metric=k, value=v) for k, v in metrics.items()]
    return rows


# storage helper funcs 
def append(path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    if path.exists():
        df = pd.concat([pd.read_parquet(path), df], ignore_index = True)
    df = harmonise(df)
    df.to_parquet(path, index = False)

def harmonise(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if df[col].dtype != object:
            continue
        types = {type(v) for v in df[col].dropna()}
        if len(types) > 1:
            df[col] = df[col].map(lambda v: v if pd.isna(v) else str(v))
    return df

def completed(path, splits = cfg.SPLITS, n_folds: int = cfg.N_FOLDS) -> set:
    df = pd.read_parquet(path)
    counts = (df[df["metric"] == "pearson_top"]
              .groupby(["spec", "model", "cell_line", "dose"])["fold"].count())
    return set(counts[counts >= len(splits) * n_folds].index)

def load_results(out_name: str = "metrics") -> pd.DataFrame:
    path = cfg.RESULTS / f"{out_name}.parquet"
    return pd.read_parquet(path) if path.exists() else pd.DataFrame()