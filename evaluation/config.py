from __future__ import annotations
from chemembed.config import DATA, ROOT
RESULTS = DATA / "results"            
FIGURES = ROOT / "figures"

RESPONSES = DATA / "responses" # Y matrices- 1 for each cell line x dose

for p in (RESPONSES, RESULTS, FIGURES):
    p.mkdir(parents = True, exist_ok = True)

##################### 
# DATASET - SCIPLEX3 
##################### 

DATASET = "sciplex3"
CELL_LINES = ("A549", "K562", "MCF7")
TIME = 24.0    
DOSES = (10.0, 100.0, 1000.0, 10000.0)
PRIMARY_DOSE = 10000.0
                       
N_HVG = 2000                          
TOP_N_DEGS = 100 # ask - should i try out a few different numbers for this?                       
MIN_CELLS_PER_CONDITION = 20
          
N_FOLDS = 5 
SPLITS = ("random", "scaffold")
CONTROL_ID_PREFIX = "CONTROL_"

SPEC_KWARGS = {
    "ecfp_raw/counts": {"with_std": False},
    "ecfp_raw/binary": {"with_std": False},
}


##################### 
# DATASET - L1000? 
##################### 

##################### 
# DATASET - Tahoe100M?
##################### 