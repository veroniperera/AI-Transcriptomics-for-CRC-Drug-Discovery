import argparse
import io

import joblib
import json
import numpy as np
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

from common import DOCK, ML, PROC, RAW, seed_all
from features import DESC_2D, DESC_3D, ad_score, bits, featurize

SORAFENIB = "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1"
CHEMBL_ID, PUBCHEM_CID = "CHEMBL1336", 216239

ap = argparse.ArgumentParser()
ap.add_argument("--chembl_sim", type=int, default=60)
ap.add_argument("--pubchem_sim", type=int, default=80)
ap.add_argument("--n_dock", type=int, default=50)
ap.add_argument("--library_csv", default=None)
args = ap.parse_args()
seed_all()


def from_chembl():
    from chembl_webresource_client.new_client import new_client
    hits = new_client.similarity.filter(smiles=SORAFENIB, similarity=args.chembl_sim)
    df = pd.DataFrame(list(hits))
    df = df[["molecule_chembl_id", "molecule_structures"]].dropna()
    df["smiles"] = df["molecule_structures"].apply(lambda d: d.get("canonical_smiles"))
    return pd.DataFrame({"id": df["molecule_chembl_id"], "smiles": df["smiles"], "source": "ChEMBL"})


def from_pubchem():
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{PUBCHEM_CID}"
           f"/property/SMILES/CSV?Threshold={args.pubchem_sim}&MaxRecords=5000")
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text))
    return pd.DataFrame({"id": "CID" + df["CID"].astype(str), "smiles": df.iloc[:, 1], "source": "PubChem"})


if args.library_csv:
    lib = pd.read_csv(args.library_csv)
else:
    lib = pd.concat([from_chembl(), from_pubchem()], ignore_index=True)
    lib.to_csv(RAW / "screening_library_raw.csv", index=False)

chooser = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()
rows = []
for r in lib.itertuples():
    m = Chem.MolFromSmiles(r.smiles) if isinstance(r.smiles, str) else None
    if m is None:
        continue
    m = uncharger.uncharge(chooser.choose(m))
    rows.append({"id": r.id, "source": r.source, "smiles": Chem.MolToSmiles(m),
                 "inchikey": inchi.MolToInchiKey(m)})
lib = pd.DataFrame(rows).drop_duplicates("inchikey").reset_index(drop=True)

params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains = FilterCatalog(params)
lib = lib[[not pains.HasMatch(Chem.MolFromSmiles(s)) for s in lib["smiles"]]].reset_index(drop=True)

train = pd.read_csv(PROC / "shmt2_curated.csv")
train_keys = set()
for s in train["smiles"]:
    m = chooser.choose(Chem.MolFromSmiles(s))
    train_keys.add(inchi.MolToInchiKey(m).split("-")[0])
lib["skeleton"] = lib["inchikey"].str.split("-").str[0]
n_overlap = int(lib["skeleton"].isin(train_keys).sum())
lib = lib[~lib["skeleton"].isin(train_keys)].reset_index(drop=True)

feat = featurize(lib["smiles"])
ok = ~feat.isna().any(axis=1)
lib, feat = lib[ok].reset_index(drop=True), feat[ok].reset_index(drop=True)

name = open(ML / "best_model.txt").read().strip()
model = joblib.load(ML / "final_model.joblib")
lib["active_prob"] = model.predict_proba(feat[DESC_2D + DESC_3D].values)[:, 1]

fp_train = np.load(PROC / "fingerprints.npz")["ecfp4"]
tr = pd.read_csv(PROC / "features.csv")
fp_tr = fp_train[(tr["split"] == "train").values]
lib["ad_score"] = ad_score(bits(lib["smiles"], "ecfp4"), fp_tr, k=5)
thr = json.load(open(ML / "applicability_domain.json"))["threshold_median_test"]
lib["in_domain"] = lib["ad_score"] >= thr

lib = lib.sort_values("active_prob", ascending=False).reset_index(drop=True)
lib.to_csv(DOCK / "screening_library_ranked.csv", index=False)
lib.head(args.n_dock).to_csv(PROC / "docking_set.csv", index=False)
print(len(lib), n_overlap, name)
