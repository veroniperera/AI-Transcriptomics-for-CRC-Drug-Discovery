import argparse

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from common import ACTIVE_PIC50, PROC, RAW, TARGET_CHEMBL

ap = argparse.ArgumentParser()
ap.add_argument("--target", default=TARGET_CHEMBL)
ap.add_argument("--raw_csv", default=None)
args = ap.parse_args()

if args.raw_csv:
    raw = pd.read_csv(args.raw_csv)
else:
    from chembl_webresource_client.new_client import new_client
    q = new_client.activity.filter(target_chembl_id=args.target, standard_type="IC50", assay_type="B")
    raw = pd.DataFrame(list(q))
    raw.to_csv(RAW / "chembl_shmt2_ic50_raw.csv", index=False)

df = raw.rename(columns={"molecule_chembl_id": "chembl_id", "canonical_smiles": "smiles"})
df = df.dropna(subset=["chembl_id", "smiles", "standard_value"])
df = df[df["standard_units"] == "nM"].copy()
df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
df = df[df["standard_value"] > 0]

mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
df = df[[m is not None for m in mols]].copy()
mols = [m for m in mols if m is not None]
df["MW"] = [Descriptors.MolWt(m) for m in mols]
df["AlogP"] = [Crippen.MolLogP(m) for m in mols]
df = df[df["MW"].between(100, 900) & df["AlogP"].between(-3, 7)]

df["pIC50"] = 9 - np.log10(df["standard_value"])
df = df.drop_duplicates(subset=["chembl_id", "smiles", "standard_value"])
df = (df.groupby(["chembl_id", "smiles"], as_index=False)
        .agg(pIC50=("pIC50", "median"), MW=("MW", "first"), AlogP=("AlogP", "first")))
df = df.drop_duplicates(subset="chembl_id", keep="first")
df["label"] = (df["pIC50"] >= ACTIVE_PIC50).astype(int)

df.to_csv(PROC / "shmt2_curated.csv", index=False)
print(len(df), df["label"].sum(), round(df["label"].mean(), 3))
