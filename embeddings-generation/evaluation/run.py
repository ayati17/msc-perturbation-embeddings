"""
Pipeline:  
    Builds (spec x regression model x cell line x split x fold) -> 1 row.
rows appended to a parquet shard as each spec finishes
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


def run(specs: list[str], models: list[str],
        cell_lines: tuple = cfg.CELL_LINES, dose: float = cfg.PRIMARY_DOSE, splits: tuple = cfg.SPLITS, 
        force: bool = False, out_name: str = "metrics") -> pd.DataFrame:

    path = cfg.RESULTS /f"{out_name}.parquet"
    done = completed(path, splits, cfg.N_FOLDS) if (path.exists() and not force) else set()

    rows = []
    for spec, cell_line in product(specs, cell_lines):
        try:
            X, Y, ids = load_xy(spec, cell_line, dose)
        except (FileNotFoundError, KeyError) as e:
            print(f"[skip] {spec} @ {cell_line}: {e}")
            continue

        print(f"\n {spec} @ {cell_line}  X{X.shape} Y{Y.shape} ")
        for model_name in models:
            if (spec, model_name, cell_line) in done:
                print(f"{model_name:10s} already done")
                continue

            started = time.perf_counter()
            new = run_one(spec, model_name, cell_line, dose, X, Y, ids, splits)
            rows.extend(new)
            score = np.mean([r["value"] for r in new if r["metric"] == "pearson_top"])
            print(f"  {model_name:10s} pearson_top = {score:+.3f}  "
                  f"({time.perf_counter() - started:.1f}s)")

        if rows:
            # adds on for each additional metric 
            append(path, rows)
            rows = []

    return load_results(out_name)


def run_one(spec, model_name, cell_line, dose, X, Y, ids, splits) -> list[dict]:
    """returns a list of rows, one per (split, fold) for the given spec/model/cell_line combo"""
    rows = []
    for split_kind in splits:
        for fold, (train, test) in enumerate(make_folds(ids, split_kind)):
            model = build(model_name, **cfg.SPEC_KWARGS.get(spec, {}))
            model.fit(X[train], Y[train])
            metrics = evaluate(Y[test], model.predict(X[test]))

            base = dict(spec = spec, model = model_name, 
                        cell_line = cell_line, dose = dose,
                        split = split_kind, fold = fold,
                        n_train = len(train), n_test = len(test), n_features = X.shape[1])
            rows += [dict(base, metric = k, value = v) for k, v in metrics.items()]
    return rows


# storage helper funcs 
def append(path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    if path.exists():
        df = pd.concat([pd.read_parquet(path), df], ignore_index = True)
    df.to_parquet(path, index = False)


def completed(path, splits = cfg.SPLITS, n_folds: int = cfg.N_FOLDS) -> set:
    """ only count a (spec, model, cell_line) combo as done if every fold is present."""
    df = pd.read_parquet(path)
    counts = (df[df["metric"] == "pearson_top"]
              .groupby(["spec", "model", "cell_line"])["fold"].count())
    return set(counts[counts >= len(splits) * n_folds].index)


def load_results(out_name: str = "metrics") -> pd.DataFrame:
    path = cfg.RESULTS / f"{out_name}.parquet"
    return pd.read_parquet(path) if path.exists() else pd.DataFrame()