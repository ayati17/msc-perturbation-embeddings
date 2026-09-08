# One module per dataset. load_raw() -> DataFrame[dataset_drug_name, smiles_raw, is_control]
from chemembed.datasets import sciplex3
# from chemembed.datasets import l1000, tahoe

LOADERS = {
    "sciplex3": sciplex3.load_raw,
    # "tahoe": tahoe.load_raw,
    # "l1000": l1000.load_raw,
}

__all__ = ["LOADERS"]