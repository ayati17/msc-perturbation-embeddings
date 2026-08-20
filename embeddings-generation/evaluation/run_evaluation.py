"""
Pipeline: 
    Choose embeddings and regression models, run the models, print results 
    
    python run_evaluation.py --check
    python run_evaluation.py --specs ecfp_raw/counts sand/z --models ridge mean #individual!
    python run_evaluation.py # runs all of SPECS x MODELS

Instructions: 
    Run with .envs/core: PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/core/bin/python -u evaluation/run_evaluation.py
"""

from __future__ import annotations
import argparse
from chemembed.artifacts import artifact_dir
from evaluation import config as cfg


# update this experiment wise: so far covered 
# ecfp + sand + string + qwen3
SPECS = [
    "ecfp_raw/counts",
    "sand/z",
    
    "string/chemberta",
    "string/selformer",
    "string/chemgpt",
    
    # to be repeated for bekko and f2llm
    "text_struct/qwen3/functional_groups",
    "text_struct/qwen3/rdkit_descriptors",
    "text_struct/qwen3/ecfp_binary",
    "text_struct/qwen3/ecfp_positional",
    "text_struct/qwen3/murcko_scaffold",
    "text_struct/qwen3/iupac",
    "text_struct/qwen3/combined",
    
    "text_semantic/qwen3/curated",
    
    "text_struct/bekko/functional_groups",
    "text_struct/bekko/rdkit_descriptors",
    "text_struct/bekko/ecfp_binary",
    "text_struct/bekko/ecfp_positional",
    "text_struct/bekko/murcko_scaffold",
    "text_struct/bekko/iupac",
    "text_struct/bekko/combined",
    
    "text_semantic/bekko/curated",
    
    "text_struct/f2llm/functional_groups",
    "text_struct/f2llm/rdkit_descriptors",
    "text_struct/f2llm/ecfp_binary",
    "text_struct/f2llm/ecfp_positional",
    "text_struct/f2llm/murcko_scaffold",
    "text_struct/f2llm/iupac",
    "text_struct/f2llm/combined",
    
    "text_semantic/f2llm/curated",
    
    
]
MODELS = ["control", "mean", "random", "ridge", "knn", "mlp", "gp"]


def main():
    args = parse_args()
    specs = args.specs or SPECS
    models = args.models or MODELS

    present = [s for s in specs if (artifact_dir(s) / "embeddings.npy").exists()]
    missing = [s for s in specs if s not in present]

    print(f"{len(present)} artifacts found, {len(missing)} missing")
    for s in missing:
        print(f"MISSING {s}")
    print(f"models: {models}")
    print(f"cell lines: {cfg.CELL_LINES} | dose: {cfg.PRIMARY_DOSE:g} nM")
    print(f"splits: {cfg.SPLITS} | folds: {cfg.N_FOLDS}")

    # x * 7 models * 3 lines * 2 splits * 5 folds 
    fits = len(present) * len(models) * len(cfg.CELL_LINES) * len(cfg.SPLITS) * cfg.N_FOLDS
    print(f"-> {fits} model fits")

    if args.check:
        return

    from evaluation.results import beats_baseline, ranking
    from evaluation.run import run

    run(present, models, force = args.force)
    print("\nranking (ridge)")
    print(ranking(model = "ridge").to_string(index = False))
    print("\nimprovement over the mean baseline")
    print(beats_baseline().to_string(index = False))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--specs", nargs ="*", help = "embedding specs; default: SPECS in here")
    p.add_argument("--models", nargs ="*", help = "model names; default: MODELS in here")
    p.add_argument("--check", action ="store_true", help = "exit")
    p.add_argument("--force", action ="store_true", help = "recompute all the rows already stored")
    return p.parse_args()

if __name__ == "__main__":
    main()