"""
Restructure the results + rankings and code for plotting annd graphs
    results: 1 row per (spec, model, cell_line, split, fold, metric)
    credit: claude for results analysis functions based on outline
    graphs, plots and funcs to add tomorrow -
        - TODO: ask how to handle dosage and eval performance - should this be added in as a variable? would need to rerun
        - model performance of MLP and ridge - tune hyperparams accordingly - TODO: what's the best way to do so? 
        - performance by cell line
        - overall highest performing model (absolute %, not compared to mean)
        - combine performance over cell line
        - compare performance of combined over individual structural embeddings
        - compare performance of the three LLMs
        - change k nearest neighbors to see if performance increases
        - compare performance of random vs scaffold
        - TODO: ask more tomorrow
        
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
    """ is the ranking stable?"""
    df = load_results() if df is None else df
    s = summary(df, metric=metric, by=("spec", "model", "cell_line"))
    return s.pivot_table(index=["spec", "model"], columns="cell_line", values="mean")


def split_gap(df: pd.DataFrame | None = None, metric: str = METRIC) -> pd.DataFrame:
    """
    Random-split score minus scaffold-split score. 
    Large gap = the embedding is leaning on scaffold similarity rather than generalising.
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
