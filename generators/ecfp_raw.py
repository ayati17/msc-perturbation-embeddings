"""
Baseline: ECFP fingerprints 
env: core
"""

from __future__ import annotations
import numpy as np
from generators._base import finish, load_inputs, parse_spec
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Hyperparameters - check if those should be changed? 
RADIUS = 2
N_BITS = 2048

def main():
    spec = parse_spec()          # variant: "counts" | "bits"
    df = load_inputs()
    
    # debugging
    # print(spec, spec.family, spec.model, spec.variant)
    # print(df.shape, df.columns.tolist())
    # print(df.head().to_string())
    
    # TODO: can implement ECFP fingerprint generation using RDKit from SMILES strings in df["smiles_canonical"] - DONE! 
    generator = rdFingerprintGenerator.GetMorganGenerator(radius = RADIUS, fpSize = N_BITS)
    smiles_list = df["smiles_canonical"].tolist()
    fingerprints = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"unparseable canonical SMILES: {smiles!r}") # should never happen </3 
        else:
            if spec.variant == "counts":
                fp = generator.GetCountFingerprintAsNumPy(mol)
            else:
                fp = generator.GetFingerprintAsNumPy(mol)
            fingerprints.append(fp)
            
    matrix = np.vstack(fingerprints).astype(np.float32)

    
    finish(spec, df["compound_id"].tolist(), matrix,
           meta = {"radius": RADIUS, "n_bits": N_BITS, "variant": spec.variant})


if __name__ == "__main__":
    main()