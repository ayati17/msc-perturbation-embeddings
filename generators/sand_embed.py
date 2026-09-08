"""
SAND - Shape-Aware Neural Descriptor (Winter et al., ICML 2026). 
env: sand

GINE encoder trained so that cosine similarity between embeddings tracks 3D shape
overlap, computed from the 2D graph from rdkit. 

  z       continuous encoder output, already L2-normalised by the model (Eqs 1-2 of
          the paper). From the no-compression checkpoint. for regression
SANDModel.load_from_checkpoints +  encode_smiles 
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from generators._base import finish, load_inputs, parse_spec

SAND_ROOT = Path(os.environ.get("SAND_ROOT", Path.home() / "SAND"))

# variant -> (checkpoint, training config, emb_key)
VARIANTS = {
    "z":      ("sand_no_comp.ckpt",   "no_compression.yml", None),
    "z_hat":  ("sand_k256_m32.ckpt",  "compression.yml",    "emb_quantized"),
    "codes":  ("sand_k256_m32.ckpt",  "compression.yml",    None),
}

BATCH_SIZE = 128


def main():
    spec = parse_spec()
    checkpoint, config, emb_key = VARIANTS[spec.variant]
    df = load_inputs()
    smiles = df["smiles_canonical"].tolist()

    model = _load_model(checkpoint, config)

    if spec.variant == "codes":
        matrix = _encode_codes(model, smiles)
    else:
        emb = model.encode_smiles(smiles, emb_key=emb_key, batch_size=BATCH_SIZE)
        matrix = emb.float().cpu().numpy()

    _report(spec, matrix)

    finish(spec, df["compound_id"].tolist(), matrix, meta={
        "checkpoint": checkpoint,
        "config": config,
        "emb_key": emb_key or "emb",
        "variant": spec.variant,
    })


def _load_model(checkpoint: str, config: str):
    import torch
    from sand.model import SANDModel

    checkpoint_path = SAND_ROOT / "checkpoints" / checkpoint
    config_path = SAND_ROOT / "sand" / "config" / config
    for p in (checkpoint_path, config_path):
        if not p.exists():
            raise FileNotFoundError(f"{p} - not found. recheck")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    map_location = None if device == "cuda" else torch.device("cpu")

    model = SANDModel.load_from_checkpoint(
        str(checkpoint_path), str(config_path), map_location=map_location
    ).eval()
    if device == "cuda":
        model = model.to("cuda")

    print(f"[sand] {checkpoint} on {device}")
    return model


def _encode_codes(model, smiles: list[str]) -> np.ndarray:
    import torch
    from sand.mol_utils.pytorch import smiles_to_batch

    out = []
    for i in range(0, len(smiles), BATCH_SIZE):
        batch = smiles_to_batch(smiles[i:i + BATCH_SIZE])
        batch = model.maybe_move_batch_to_gpu(batch)
        with torch.no_grad():
            batch = model.backbone.encode(batch)
            batch = model.backbone.forward_batch_ivfpq(batch)

        code_keys = [k for k, v in batch.items()
                     if torch.is_tensor(v) and not v.is_floating_point() and v.dim() <= 2]
        if not code_keys:
            raise KeyError(
                f"no integer code tensors in batch; keys were {sorted(batch.keys())}. "
                "Inspect forward_batch_ivfpq and pick the IVF/PQ index tensors by hand."
            )
        print(f"[sand] code keys: {code_keys}")
        codes = torch.cat([batch[k].reshape(len(smiles[i:i + BATCH_SIZE]), -1)
                           for k in sorted(code_keys)], dim=1)
        out.append(codes.cpu().numpy().astype(np.int32))

    return np.vstack(out)


def _report(spec, matrix: np.ndarray) -> None:
    print(f"[{spec}] {matrix.shape} {matrix.dtype}")
    if matrix.dtype.kind == "f":
        norms = np.linalg.norm(matrix, axis=1)
        print(f"[{spec}] norms: min {norms.min():.4f} max {norms.max():.4f} "
              f"(z is L2-normalised by the model, so expect ~1.0)")


if __name__ == "__main__":
    main()