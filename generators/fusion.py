"""
Fusion of multiple embedding blocks into one artifact. 
Reads in existing artifacts by compound_id and writes out a new artificat for a fusion embedding 
env: core

Run after all the individual embeddings have been evaluated 

Selected inputs: 
    BASE_STRING - best of the SMILES / SELFIES sequence models
    TEXT_MODEL - chosen text embedders (bekko / f2llm / qwen3)


FUSION COMBINATIONS (check google doc)
Two bases:                          A = BASE_STRING              (sequence model)
                                    B = sand/z                   (3D shape)

Seven individual partners:          text_struct/<M>/functional_groups
                                    text_struct/<M>/rdkit_descriptors
                                    text_struct/<M>/ecfp_binary
                                    text_struct/<M>/ecfp_positional
                                    text_struct/<M>/murcko_scaffold
                                    text_struct/<M>/iupac
                                    text_semantic/<M>/curated

  2 bases x 7 partners                                          = 14   (2 blocks each)

Combined structural text (note that the text was joined as TEXT, not as embeddings)
    A x text_struct/<M>/combined                                =  1   (2 blocks)
    B x text_struct/<M>/combined                                =  1   (2 blocks)

Combined structural text plus semantic:
    A x text_struct/<M>/combined x text_semantic/<M>/curated    =  1   (3 blocks)
    B x text_struct/<M>/combined x text_semantic/<M>/curated    =  1   (3 blocks)

                                                          total = 18 possible combos
                                                        

Spec grammar:  fusion/<method>/<combo>      e.g.  fusion/pca_concat/sand_x_iupac
    pca_concat  standardise + PCA each block to N_COMPONENTS first, so blocks enter
                on comparable footing. Likely the strongest practical fusion here.
    additive kernel GP - complete this one out tomorrow 
"""

from __future__ import annotations

import numpy as np

from generators._base import finish, load_inputs, parse_spec
from chemembed.registry import COMBOS
from chemembed.registry import FUSION_METHODS

# N_COMPONENTS = 50                     
def main():
    spec = parse_spec()                       
    method, combo = spec.model, spec.variant
    blocks = COMBOS[combo]

    df = load_inputs()
    ids = df["compound_id"].tolist()

    matrices = [_load_block(b, ids) for b in blocks]
    print(f"[{spec}] {len(blocks)} blocks: "
          f"{[f'{b} {m.shape[1]}d' for b, m in zip(blocks, matrices)]}")

    matrix = FUSERS[method](matrices)

    finish(spec, ids, matrix, meta={
        "method": method,
        "combo": combo,
        "blocks": list(blocks),
        "block_dims": [int(m.shape[1]) for m in matrices],
        "n_components": N_COMPONENTS if method != "concat" else None,
    })


def _load_block(spec: str, ids: list[str]) -> np.ndarray:
    # note: RAW
    from chemembed.artifacts import load_embedding
    return load_embedding(spec, ids, normalize=False, standardize=False)



def fuse_concat(matrices: list[np.ndarray]) -> np.ndarray:
    return np.hstack(matrices).astype(np.float32)

from functools import partial

def fuse_pca_concat(matrices: list[np.ndarray], k: int) -> np.ndarray:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    reduced = []
    for block in matrices:
        n_components = min(k, block.shape[1], len(block) - 1)
        pca = PCA(n_components=n_components, random_state=17)
        reduced.append(pca.fit_transform(StandardScaler().fit_transform(block)))
        print(f"[fusion]   {block.shape[1]:5d}d ---> {n_components:3d} comps "
              f"({pca.explained_variance_ratio_.sum():.1%} variance)")
    return np.hstack(reduced).astype(np.float32)


FUSERS = {
    # "concat": fuse_concat,
    # "pca_concat100": partial(fuse_pca_concat, k = 100),
    "pca_concat50": partial(fuse_pca_concat, k = 50),
}
assert set(FUSERS) == set(FUSION_METHODS), (
    f"registry says {FUSION_METHODS}, fusion.py has {tuple(FUSERS)}")


if __name__ == "__main__":
    main()