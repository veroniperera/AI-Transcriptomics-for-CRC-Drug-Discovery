import numpy as np
import pandas as pd
from rdkit.Chem.Scaffolds import MurckoScaffold

from common import PROC, SEED, seed_all
from features import bits, featurize

seed_all()
df = pd.read_csv(PROC / "shmt2_curated.csv")

df["scaffold"] = [MurckoScaffold.MurckoScaffoldSmiles(smiles=s) for s in df["smiles"]]
groups = df.groupby("scaffold").indices
scaffolds = sorted(groups)
rng = np.random.RandomState(SEED)
perm = rng.permutation(len(scaffolds))
n_test = int(round(0.2 * len(df)))
test_idx = []
for i in perm:
    if len(test_idx) >= n_test:
        break
    test_idx.extend(groups[scaffolds[i]])
df["split"] = "train"
df.loc[df.index[test_idx], "split"] = "test"

feat = featurize(df["smiles"])
out = pd.concat([df.reset_index(drop=True), feat], axis=1)
out.to_csv(PROC / "features.csv", index=False)

np.savez_compressed(PROC / "fingerprints.npz",
                    ecfp4=bits(df["smiles"], "ecfp4"),
                    ecfp6=bits(df["smiles"], "ecfp6"),
                    rdkit=bits(df["smiles"], "rdkit"))

print(len(scaffolds), (out["split"] == "train").sum(), (out["split"] == "test").sum(),
      set(out.loc[out.split == "train", "scaffold"]) & set(out.loc[out.split == "test", "scaffold"]),
      int(feat.isna().sum().sum()))
