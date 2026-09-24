import os
import random
import warnings
from pathlib import Path

import numpy as np

SEED = 42
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
ML = ROOT / "results" / "ml"
DOCK = ROOT / "results" / "docking"
DUDE = ROOT / "results" / "dude"
for p in (RAW, PROC, ML, DOCK, DUDE):
    p.mkdir(parents=True, exist_ok=True)

TARGET_CHEMBL = "CHEMBL4295747"
ACTIVE_PIC50 = 5.35
PDB_ID = "7BYI"
BOX_CENTER = (9.55, -69.71, 1.65)
BOX_SIZE = (22.5, 22.5, 36.6)


def seed_all(seed=SEED):
    warnings.filterwarnings("ignore")
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
