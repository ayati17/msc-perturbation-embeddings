"""Driver for the CellFlow backbone evaluation.

    python run_backbone.py --check
    python run_backbone.py --specs ecfp_raw/counts text_semantic/qwen3/curated --folds 0
    python run_backbone.py --iterations 2000 --folds 0      #  test
    
    CUDA_VISIBLE_DEVICES=3 PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/cellflow/bin/python -u backbone_eval/run_backbone.py \
    --specs sand/z \
    --folds 0 1 2 3 4 \
    --iterations 20_000

Run with the env that has cellflow + jax installed:

    PYTHONPATH=. CHEMEMBED_ROOT=$(pwd) .envs/cellflow/bin/python -u run_backbone.py
"""

from __future__ import annotations

import argparse

from backbone_eval import config as cfg
from backbone_eval.reps import available

# Selected embeddings from individual and fusion
# Each entry is a full CellFlow training run per fold.
SPECS = [
    "ecfp_raw/counts",                  # baseline - CellFlow's own default input
    
    # best embeddings from the regression models 
    # to be filled :DD 

]


def main():
    args = _parse_args()
    specs = available(args.specs or SPECS)
    folds = tuple(args.folds) if args.folds else tuple(range(cfg.N_FOLDS))

    print(f"specs:      {specs}")
    print(f"folds:      {folds}")
    print(f"iterations: {args.iterations}")
    print(f"-> {len(specs) * len(folds)} CellFlow training runs")

    if args.check:
        return

    from backbone_eval.run import run, summary

    run(specs, folds=folds, n_iterations=args.iterations, force=args.force)
    print("\n" + summary().to_string(index=False))


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--specs", nargs="*")
    p.add_argument("--folds", nargs="*", type=int)
    p.add_argument("--iterations", type=int, default=cfg.N_ITERATIONS)
    p.add_argument("--check", action="store_true")
    p.add_argument("--force", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    main()