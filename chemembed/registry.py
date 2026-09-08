"""
kinds of embeddings that can be generated - by which script, in which env
setup: <family>[/<model-name>][/<variant>]

Samples!
    ecfp_raw/counts
    sand/z
    string/chemberta
    text_struct/qwen3/iupac
    text_semantic/bekko
"""

from __future__ import annotations
import fnmatch
from dataclasses import dataclass, field


@dataclass(frozen = True)
class TextModel:
    hf_id: str
    pooling: str = "mean"           # "mean" for encoders, "last" for decoder embedders
    padding_side: str = "right"     # must be "left" when pooling == "last"
    max_length: int = 8192
    prefix: str = ""                


# Text embedding models, shared by the text_structural and text_semantic setups
TEXT_MODELS: dict[str, TextModel] = {
    # DISCARD LLAMA
    "llama": TextModel("meta-llama/Llama-3.1-8B", pooling = "mean", max_length = 8192), #TODO - hf token kinda buggy
    "bekko": TextModel("hotchpotch/bekko-embedding-v1-a25m", pooling = "mean", max_length = 8192),
    "f2llm": TextModel("codefuse-ai/F2LLM-v2-0.6B", pooling = "last", padding_side = "left"),
    "qwen3": TextModel(
        "Qwen/Qwen3-Embedding-0.6B",
        pooling = "last",
        padding_side = "left",
        prefix = "Instruct: Represent this molecule description for retrieval\nQuery: ",
    ),
    # "f2llm-14b": TextModel("codefuse-ai/F2LLM-v2-14B", pooling = "last", padding_side = "left"),
}


@dataclass(frozen=True)
class StringModel:
    hf_id: str | None
    input: str # "smiles" or "selfies"
    script: str
    env: str
    pooling: str = "mean"
    padding_side: str = "right"     # has gotta be "left" when pooling == "last"
    max_length: int = 512           # check if i'd need to extend this 
    prefix: str = ""


STRING_MODELS: dict[str, StringModel] = {
    # SMILES Models
    "chemberta": StringModel("DeepChem/ChemBERTa-77M-MTR", "smiles", "generators/string_hf.py", "hf"),
    # "simson": StringModel(None, "smiles", "generators/string_simson.py", "simson"), # uhh tbd 
    
    # SELFIES models
    "selformer": StringModel("HUBioDataLab/SELFormer", "selfies","generators/string_hf.py", "hf"),
    "chemgpt": StringModel("ncfrey/ChemGPT-1.2B", "selfies", "generators/string_hf.py", 
                           "hf", pooling = "last", padding_side = "left"),
    
}

# Structural Fields 
STRUCT_FIELDS = (
    "functional_groups",
    "rdkit_descriptors",
    "ecfp_binary",
    "ecfp_positional",
    "murcko_scaffold",
    "iupac",
    "combined", # all 6 combined as the input text
)

# adding in for the fusion embeddings 
# based on the individual embeddings sweep 
BASE_STRING = "string/chemberta"      # chemberta (SMILES) > selformer (SELFIES)
TEXT_MODEL = "f2llm"                  # bekko / f2llm > qwen3


BASES = {
    "str": BASE_STRING,
    "sand": "sand/z",
}

PARTNERS = {
    "fgroups":  f"text_struct/{TEXT_MODEL}/functional_groups",
    "rdkit":    f"text_struct/{TEXT_MODEL}/rdkit_descriptors",
    "ecfpbin":  f"text_struct/{TEXT_MODEL}/ecfp_binary",
    "ecfppos":  f"text_struct/{TEXT_MODEL}/ecfp_positional",
    "murcko":   f"text_struct/{TEXT_MODEL}/murcko_scaffold",
    "iupac":    f"text_struct/{TEXT_MODEL}/iupac",
    "semantic": f"text_semantic/{TEXT_MODEL}/curated",
}

COMBINED_STRUCT = f"text_struct/{TEXT_MODEL}/combined"
SEMANTIC = f"text_semantic/{TEXT_MODEL}/curated"


def _build_combos() -> dict[str, tuple[str, ...]]:
    combos: dict[str, tuple[str, ...]] = {}

    # 14: each base with each individual partner
    for base_name, base_spec in BASES.items():
        for partner_name, partner_spec in PARTNERS.items():
            combos[f"{base_name}_x_{partner_name}"] = (base_spec, partner_spec)

    # 2: each base with the combined structural text
    for base_name, base_spec in BASES.items():
        combos[f"{base_name}_x_combined"] = (base_spec, COMBINED_STRUCT)

    # 2: each base with combined structural text and semantic (three blocks)
    for base_name, base_spec in BASES.items():
        combos[f"{base_name}_x_combined_semantic"] = (base_spec, COMBINED_STRUCT, SEMANTIC)

    return combos


COMBOS = _build_combos()
FUSION_METHODS = ("concat", "pca_concat50", "pca_concat100")


@dataclass(frozen = True)
# Look at Google Doc for overall organization
# TODO: ask Jonathan for additional baseline if needed 
class Family:
    script: str
    env: str
    models: tuple[str, ...] = ()        
    variants: tuple[str, ...] = ()      
    per_model: dict = field(default_factory=dict)  # model -> (script, env) override


FAMILIES: dict[str, Family] = {
    # Baseline 
    "ecfp_raw": Family(                                 
        "generators/ecfp_raw.py", 
        "core",
        variants = ("counts", "bits"),
    ),
    "sand": Family(
        "generators/sand_embed.py", "sand",
        variants = ("z", "z_hat"),
    ),
    "string": Family(
        "generators/string_hf.py", 
        "hf",
        models = tuple(STRING_MODELS),
        per_model = {k: (m.script, m.env) for k, m in STRING_MODELS.items()},
    ),
    "text_struct": Family(
        "generators/text_struct.py", 
        "hf",
        models = tuple(TEXT_MODELS),
        variants = STRUCT_FIELDS,
    ),
    "text_semantic": Family(
        "generators/text_semantic.py", "hf",
        models = tuple(TEXT_MODELS),
        variants = ("chembl", "curated"),
    ),
    # "graph": Family(
    #     "generators/graph_preprct.py", 
    #     "graph",
    #     models = ("preprct",),
    # ),
    "fusion": Family(
        "generators/fusion.py", "core",
        models=FUSION_METHODS, #TODO: add to this
        variants=tuple(COMBOS),
    ),
    
}


@dataclass(frozen = True)
class Spec:
    family: str
    model: str | None = None
    variant: str | None = None

    @classmethod
    def parse(cls, s: str) -> "Spec":
        parts = s.strip("/").split("/")
        family, rest = parts[0], parts[1:]
        if family not in FAMILIES:
            raise ValueError(f"unknown family {family!r} in spec {s!r}")
        fam = FAMILIES[family]
        model = variant = None
        for part in rest:
            if part in fam.models and model is None:
                model = part
            elif part in fam.variants and variant is None:
                variant = part
            else:
                raise ValueError(f"cannot place {part!r} in spec {s!r}")
        if fam.models and model is None:
            raise ValueError(f"spec {s!r} needs a model: one of {fam.models}")
        if fam.variants and variant is None:
            raise ValueError(f"spec {s!r} needs a variant: one of {fam.variants}")
        return cls(family, model, variant)

    def __str__(self) -> str:
        return "/".join(p for p in (self.family, self.model, self.variant) if p)

    @property
    def script(self) -> str:
        fam = FAMILIES[self.family]
        return fam.per_model.get(self.model, (fam.script, fam.env))[0]

    @property
    def env(self) -> str:
        fam = FAMILIES[self.family]
        return fam.per_model.get(self.model, (fam.script, fam.env))[1]


def all_specs() -> list[Spec]:
    out = []
    for family, fam in FAMILIES.items():
        for model in fam.models or (None,):
            for variant in fam.variants or (None,):
                out.append(Spec(family, model, variant))
    return out


def expand_specs(patterns: list[str]) -> list[Spec]:
    """Expand glob patterns against the full spec list, preserving order."""
    seen, out = set(), []
    universe = all_specs()
    for pattern in patterns:
        matched = [s for s in universe if fnmatch.fnmatch(str(s), pattern)]
        if not matched:
            matched = [Spec.parse(pattern)]     
        for spec in matched:
            if str(spec) not in seen:
                seen.add(str(spec))
                out.append(spec)
    return out