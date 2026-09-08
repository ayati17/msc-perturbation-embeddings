
"""
SMILES -> the six structural text renderings, cached to JSON. 
 functional groups (SMARTS catalogue → prose), 
 RDKit descriptors (selected values → labelled text), 
 ECFP binary (on-bits as tokens), 
 ECFP positional (on-bits + atom environments), 
 Murcko scaffold, 
 IUPAC (pulled from the PubChem cache). 
 combined = concatenation of the other six, in FIELD_ORDER.
env: core
Output: data cached to data/cache/descriptions.json    {compound_id: {field: text}}
Steps: 
1. Run once after build_compounds() and build_pubchem(), BEFORE any text spec:
    build_descriptions()
"""
from __future__ import annotations
import json
import numpy as np
from chemembed.config import CACHE
CACHE_FILE = CACHE / "descriptions.json"

FIELD_ORDER = (
    "functional_groups",
    "rdkit_descriptors",
    "ecfp_binary",
    "ecfp_positional",
    "murcko_scaffold",
    "iupac",
)

FUNC_GROUPS = {
    "fr_benzene": "benzene ring",
    "fr_NH0": "tertiary amine",
    "fr_NH1": "secondary amine",
    "fr_NH2": "primary amine",
    "fr_C_O": "carbonyl group",
    "fr_Ar_N": "aromatic nitrogen",
    "fr_amide": "amide",
    "fr_priamide": "primary amide",
    "fr_bicyclic": "bicyclic ring system",
    "fr_aniline": "aniline",
    "fr_Ar_NH": "aromatic secondary amine",
    "fr_ether": "ether",
    "fr_methoxy": "methoxy group",
    "fr_halogen": "halogen",
    "fr_alkyl_halide": "alkyl halide",
    "fr_aryl_methyl": "aryl methyl group",
    "fr_pyridine": "pyridine ring",
    "fr_imidazole": "imidazole ring",
    "fr_Nhpyrrole": "pyrrole NH",
    "fr_piperdine": "piperidine ring",
    "fr_piperzine": "piperazine ring",
    "fr_morpholine": "morpholine ring",
    "fr_phenol": "phenol",
    "fr_Al_OH": "aliphatic hydroxyl",
    "fr_ketone": "ketone",
    "fr_nitrile": "nitrile",
    "fr_COO": "carboxylic acid",
    "fr_urea": "urea",
    "fr_sulfonamd": "sulfonamide",
    "fr_sulfone": "sulfone",
    "fr_unbrch_alkane": "unbranched alkane chain",
    "fr_N_O": "N-hydroxy group (hydroxamic acid or hydroxylamine)",
}

RADIUS = 2
N_BITS = 2048


def build_descriptions(force: bool = False) -> dict:
    from chemembed.compounds import load_all
    from chemembed.text.pubchem import load_iupac_names

    frags = _load_fragments()
    iupac = load_iupac_names()

    all_compounds = load_all()
    keep = ~all_compounds["is_control"].astype(bool) & all_compounds["smiles_canonical"].notna()
    df = all_compounds[keep]

    out = {}
    for row in df.itertuples(index = False):
        smiles = row.smiles_canonical
        fields = {
            "functional_groups": _functional_groups(smiles, frags),
            "rdkit_descriptors": _rdkit_descriptors(smiles),
            "ecfp_binary": _ecfp_binary(smiles),
            "ecfp_positional": _ecfp_positional(smiles),
            "murcko_scaffold": _murcko_scaffold(smiles),
            "iupac": iupac.get(row.compound_id, ""),
        }
        fields["combined"] = "\n".join(fields[f] for f in FIELD_ORDER if fields[f])
        out[row.compound_id] = fields

    CACHE_FILE.write_text(json.dumps(out, indent=2))
    _report(out)
    return out


def load_descriptions() -> dict:
    """called by generators"""
    return json.loads(CACHE_FILE.read_text())


def _load_fragments() -> dict:
    from rdkit.Chem import Fragments
    frags = {name: fn for name, fn in Fragments.__dict__.items()
             if name.startswith("fr_") and callable(fn)}
    missing = [f for f in FUNC_GROUPS if f not in frags]
    assert not missing, f"Missing RDKit fragments: {missing}"
    return frags


# Each helper func returns a string for the LLMs to embed.

def _functional_groups(smiles: str, frags: dict) -> str:
    """'2 amides, 1 nitrile, 3 benzene rings'"""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""

    present = [
        f"{n} {_plural(label, n)}"
        for fragment, label in FUNC_GROUPS.items()
        if (n := frags[fragment](mol)) > 0
    ]
    return ", ".join(present) if present else "No functional groups detected."


def _rdkit_descriptors(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""

    descriptors = {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "num_h_donors": rdMolDescriptors.CalcNumHBD(mol),
        "num_h_acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "num_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "heavy_atom_count": rdMolDescriptors.CalcNumHeavyAtoms(mol),
        "num_stereocentres": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "qed": Descriptors.qed(mol),
    }
    return "; ".join(
        f"{k} {v:.2f}" if isinstance(v, float) else f"{k} {v}"
        for k, v in descriptors.items()
    )


def _ecfp_binary(smiles: str) -> str:
    """ Returns a text string like: 'ECFP4 binary: 100010...' """
    arr = _ecfp_array(smiles)
    if arr is None:
        return ""
    return "ECFP4 binary: " + "".join(arr.astype(str).tolist())


def _ecfp_positional(smiles: str) -> str:
    """ Returns a text string like: 'ECFP4 active bit positions: 0, 5, 7, 9'
    """
    arr = _ecfp_array(smiles)
    if arr is None:
        return ""
    on_bits = np.flatnonzero(arr)
    if not len(on_bits):
        return "ECFP4 active bit positions: none"
    return "ECFP4 active bit positions: " + ", ".join(map(str, on_bits.tolist()))


def _ecfp_array(smiles: str) -> np.ndarray | None:
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius = RADIUS, fpSize = N_BITS)
    return generator.GetFingerprintAsNumPy(mol)


def _murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        return f"Murcko scaffold: {scaffold}" if scaffold else "Murcko scaffold: acyclic"
    except Exception:
        return ""


def _plural(label: str, n: int) -> str:
    return label if n == 1 else label + "s"


def _report(out: dict) -> None:
    print(f"[descriptions] {len(out)} compounds -> {CACHE_FILE}")
    for field in FIELD_ORDER + ("combined",):
        empty = sum(1 for v in out.values() if not v[field].strip())
        lengths = sorted(len(v[field]) for v in out.values())
        print(f"{field:20s} empty={empty:3d}  median={lengths[len(lengths) // 2]:5d} chars")