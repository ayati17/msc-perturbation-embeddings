""" CellFlow backbone evaluation configurations """

from __future__ import annotations
from chemembed.config import DATA, ROOT

RESULTS = DATA / "results" / "backbone"
FIGURES = ROOT / "figures"
for _p in (RESULTS, FIGURES):
    _p.mkdir(parents = True, exist_ok = True)


# Dataset related 
DATASET = "sciplex3"
TIME = 24.0
N_HVG = 2000
CONTROL_TOKEN = "control"          

N_FOLDS = 5
SEED = 17


# Cellflow specific 
N_PCA = 100                        
N_ITERATIONS = 20000
BATCH_SIZE = 256

# reused for every embedding 
# TODO: check the cellflow documentation for the layers before pool and layers after pool, 
# one of them's slowing me down - gotta delete 
MODEL_KWARGS = dict(
    condition_mode = "deterministic",
    regularization = 0.0,
    pooling = "mean", #original - attention_token
    layers_before_pool = {
        "drug": {"layer_type": "mlp", "dims": [256, 256], "dropout_rate": 0.0},
        "log_dose": {"layer_type": "mlp", "dims": [32, 32], "dropout_rate": 0.0},
    },
    layers_after_pool = { "layer_type": "mlp", "dims": [128], "dropout_rate": 0.0}, # [256, 256]
    condition_embedding_dim = 128,
    cond_output_dropout = 0.2,
    pool_sample_covariates = False,
    probability_path = {"constant_noise": 0.5},
)

BASELINE = "ecfp_raw/counts" 
# CellFlow's default is Morgan fingerprints, todo - justify equivalence in the paper