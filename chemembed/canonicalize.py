"""
Structure standardisation. 

Steps 
1. Run once, per dataset 
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_UNCHARGER = rdMolStandardize.Uncharger()
_CHOOSER = rdMolStandardize.LargestFragmentChooser()


@dataclass(frozen = True)
class Standardized:
    smiles: str
    inchikey: str


def standardize(smiles: str) -> Standardized | None:
    """
    If RDKIT can't parse the SMILES string --> None; 
    callers must handle it.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = _CHOOSER.choose(mol)          # strip salts / counterions
    mol = _UNCHARGER.uncharge(mol)
    mol = rdMolStandardize.Normalize(mol)
    Chem.AssignStereochemistry(mol, cleanIt = True, force = True)
    canonical = Chem.MolToSmiles(mol, isomericSmiles = True)
    key = Chem.MolToInchiKey(mol)
    if not key:
        return None
    return Standardized(canonical, key)