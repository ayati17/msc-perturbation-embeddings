"""
Restructure the results + rankings and code for plotting annd graphs
    results: 1 row per (spec, model, cell_line, split, fold, metric)
    credit: claude for results analysis functions based on outline
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from evaluation import config as cfg
from evaluation.run import load_results

METRIC = "pearson_top"          


def summary(df: pd.DataFrame | None = None, metric: str = METRIC,
            by: tuple = ("spec", "model")) -> pd.DataFrame:
    """ 
    mean, sd and a 95% interval over folds/cell lines.
    """
    df = load_results() if df is None else df
    sub = df[df["metric"] == metric]
    out = (sub.groupby(list(by))["value"]
              .agg(mean="mean", sd="std", n="count")
              .reset_index())
    out["ci95"] = 1.96 * out["sd"] / np.sqrt(out["n"])
    return out.sort_values("mean", ascending=False)


def ranking(df: pd.DataFrame | None = None, model: str = "ridge",
            metric: str = METRIC) -> pd.DataFrame:
    """ embeddings ranked for one model, with the baselines kept in for reference."""
    df = load_results() if df is None else df
    keep = df["model"].isin([model, "mean", "control", "random"])
    return summary(df[keep], metric=metric, by=("spec", "model"))


def beats_baseline(df: pd.DataFrame | None = None, metric: str = METRIC) -> pd.DataFrame:
    """
    checks if embedding beats predicting the mean 
    compares (spec, model) to mean baseline fold-for-fold
    """
    df = load_results() if df is None else df
    sub = df[df["metric"] == metric]

    keys = ["cell_line", "split", "fold"]
    baseline = (sub[sub["model"] == "mean"]
                .groupby(keys)["value"].mean().rename("baseline"))

    merged = sub[sub["model"] != "mean"].merge(baseline, on=keys, how="left")
    merged["delta"] = merged["value"] - merged["baseline"]

    out = (merged.groupby(["spec", "model"])["delta"]
                 .agg(mean="mean", sd="std", n="count")
                 .reset_index())
    out["ci95"] = 1.96 * out["sd"] / np.sqrt(out["n"])
    out["beats_mean"] = out["mean"] - out["ci95"] > 0
    return out.sort_values("mean", ascending=False)


def by_cell_line(df: pd.DataFrame | None = None, metric: str = METRIC) -> pd.DataFrame:
    # todo: FURTHER ANALYSIS for cell line generalizability? is the ranking stable?
    df = load_results() if df is None else df
    s = summary(df, metric=metric, by=("spec", "model", "cell_line"))
    return s.pivot_table(index=["spec", "model"], columns="cell_line", values="mean")


def split_gap(df: pd.DataFrame | None = None, metric: str = METRIC) -> pd.DataFrame:
    """
    Random-split score minus scaffold-split score. 
    """
    df = load_results() if df is None else df
    s = summary(df, metric=metric, by=("spec", "model", "split"))
    wide = s.pivot_table(index=["spec", "model"], columns="split", values="mean")
    if {"random", "scaffold"} <= set(wide.columns):
        wide["gap"] = wide["random"] - wide["scaffold"]
    return wide.sort_values("gap", ascending=False) if "gap" in wide else wide


# plots and figures
def plot_ranking(df: pd.DataFrame | None = None, model: str = "ridge",
                 metric: str = METRIC, save: bool = True):
    import matplotlib.pyplot as plt

    tab = ranking(df, model=model, metric=metric).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(tab) + 1.5))
    colours = ["0.6" if m in ("mean", "control", "random") else "C0"
               for m in tab["model"]]
    ax.barh(tab["spec"] + " / " + tab["model"], tab["mean"],
            xerr=tab["ci95"], color=colours)
    ax.set_xlabel(metric)
    ax.set_title(f"{metric} by embedding ({model}); grey = baselines")
    fig.tight_layout()
    if save:
        path = cfg.FIGURES / f"ranking_{model}_{metric}.png"
        fig.savefig(path, dpi=150)
        print(f"saved {path}")
    return fig

def headline(df=None, metric="pearson_top"):
    df = load_results() if df is None else df
    abs_ = summary(df, metric=metric).rename(columns={"mean": "absolute"})
    rel = beats_baseline(df, metric=metric).rename(columns={"mean": "vs_mean"})
    out = abs_.merge(rel[["spec", "model", "vs_mean", "ci95", "beats_mean"]],
                     on=["spec", "model"], suffixes=("_abs", "_rel"))
    return out.sort_values("vs_mean", ascending=False)

def all_metrics(df=None, models: tuple | None = None) -> pd.DataFrame:
    """Wide table: rows are (spec, model), columns are every metric."""
    df = load_results() if df is None else df
    if models:
        df = df[df["model"].isin(models)]
    return (df.groupby(["spec", "model", "metric"])["value"]
              .mean().unstack("metric").round(3))


def report(df=None, metric: str = METRIC, save: bool = True) -> str:
    df = load_results() if df is None else df
    if df.empty:
        return "no results yet"

    parts = [
        f"# Results report ({metric})",
        f"\nspecs: {df['spec'].nunique()} | models: {df['model'].nunique()} | "
        f"rows: {len(df)}",
        "\n\n## Headline - absolute score and improvement over the mean baseline",
        headline(df, metric).to_string(index=False),
        "\n\n## All metrics (mean over folds and cell lines)",
        all_metrics(df).to_string(),
        "\n\n## By cell line - is the ranking stable?",
        by_cell_line(df, metric).round(3).to_string(),
        "\n\n## Random vs scaffold split - large gap means scaffold memorisation",
        split_gap(df, metric).round(3).to_string(),
    ]
    text = "\n".join(parts)
    print(text)

    if save:
        path = cfg.RESULTS / f"report_{metric}.md"
        path.write_text(text)
        print(f"\nsaved {path}")
    return text


def hyperparams(df=None, metric: str = METRIC, params: tuple = (),
                by_spec: bool = True) -> pd.DataFrame:
    """
    One row per hyperparameter setting: test score, train score, and the gap.
    """
    df = load_results() if df is None else df
    if not params:
        params = tuple(c for c in df.columns
                       if c.startswith("p_") and df[c].nunique() > 1)

    keys = (["spec"] if by_spec else []) + ["model", *params]
    test = (df[df["metric"] == metric].groupby(keys, dropna=False)["value"]
              .agg(test="mean", sd="std", n="count").reset_index())
    train = (df[df["metric"] == f"train_{metric}"].groupby(keys, dropna=False)["value"]
               .mean().rename("train").reset_index())

    out = test.merge(train, on=keys, how="left")
    out["gap"] = out["train"] - out["test"]
    out["ci95"] = 1.96 * out["sd"] / np.sqrt(out["n"])
    return out.sort_values("test", ascending=False)


def appendix_table(out_name: str = "mlp_sweep", metric: str = METRIC) -> pd.DataFrame:
    """
    Every hyperparameter setting ever run, for the supplementary info.
    """
    from evaluation.run import load_results
    df = load_results(out_name)
    params = [c for c in df.columns if c.startswith("p_") and df[c].nunique() > 1]

    test = (df[df["metric"] == metric].groupby(params)["value"]
              .agg(test="mean", sd="std", n="count").reset_index())
    train = (df[df["metric"] == f"train_{metric}"].groupby(params)["value"]
               .mean().rename("train").reset_index())

    out = test.merge(train, on=params, how="left")
    out["gap"] = out["train"] - out["test"]
    out["ci95"] = 1.96 * out["sd"] / np.sqrt(out["n"])
    return out.sort_values("test", ascending=False).round(4)
