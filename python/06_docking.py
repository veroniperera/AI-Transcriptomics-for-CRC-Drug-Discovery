import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolAlign, rdMolDescriptors

from common import BOX_CENTER, BOX_SIZE, DOCK, PDB_ID, PROC, RAW, SEED, seed_all
from features import embed

ap = argparse.ArgumentParser()
ap.add_argument("--library", default=str(PROC / "docking_set.csv"))
ap.add_argument("--cpu", type=int, default=8)
ap.add_argument("--redock", default=None)
args = ap.parse_args()
seed_all()

WORK = DOCK / "work"
WORK.mkdir(exist_ok=True)
BOX = ["--center_x", str(BOX_CENTER[0]), "--center_y", str(BOX_CENTER[1]), "--center_z", str(BOX_CENTER[2]),
       "--size_x", str(BOX_SIZE[0]), "--size_y", str(BOX_SIZE[1]), "--size_z", str(BOX_SIZE[2])]


def prepare_receptor():
    pdb = RAW / f"{PDB_ID}.pdb"
    if not pdb.exists():
        pdb.write_text(requests.get(f"https://files.rcsb.org/download/{PDB_ID}.pdb", timeout=120).text)
    prot = WORK / "receptor.pdb"
    prot.write_text("".join(l for l in pdb.read_text().splitlines(True) if l.startswith("ATOM")) + "END\n")
    out = WORK / "receptor.pdbqt"
    subprocess.run(["obabel", str(prot), "-O", str(out), "-xr", "-h", "--partialcharge", "gasteiger"], check=True)
    return out


def props(mol):
    return {k: mol.GetProp(k) for k in mol.GetPropNames()}


def smina_dock(receptor, lig_sdf, out_sdf):
    subprocess.run(["smina", "-r", str(receptor), "-l", str(lig_sdf), "-o", str(out_sdf), *BOX,
                    "--exhaustiveness", "16", "--num_modes", "9", "--seed", str(SEED),
                    "--cpu", str(args.cpu)], check=True, capture_output=True)
    return [m for m in Chem.SDMolSupplier(str(out_sdf)) if m is not None]


def gnina_rescore(receptor, pose_sdf, out_sdf):
    subprocess.run(["gnina", "-r", str(receptor), "-l", str(pose_sdf), "--score_only",
                    "--cnn_scoring", "rescore", "-o", str(out_sdf), "--seed", str(SEED),
                    "--cpu", str(args.cpu)], check=True, capture_output=True)
    return next(m for m in Chem.SDMolSupplier(str(out_sdf)) if m is not None)


receptor = prepare_receptor()

if args.redock:
    ref = Chem.SDMolSupplier(args.redock, removeHs=True)[0]
    poses = smina_dock(receptor, args.redock, WORK / "redock_out.sdf")
    rmsd = rdMolAlign.CalcRMS(Chem.RemoveHs(poses[0]), ref)
    (DOCK / "redocking_rmsd.txt").write_text(f"{rmsd:.2f}\n")
    print("re-docking RMSD:", round(rmsd, 2))

lib = pd.read_csv(args.library)
rows = []
for r in lib.itertuples():
    m = Chem.MolFromSmiles(r.smiles)
    mh = embed(m)
    if mh is None:
        continue
    mh.SetProp("_Name", str(r.id))
    inp, out, resc = WORK / f"{r.id}.sdf", WORK / f"{r.id}_smina.sdf", WORK / f"{r.id}_gnina.sdf"
    w = Chem.SDWriter(str(inp))
    w.write(mh)
    w.close()
    try:
        best = smina_dock(receptor, inp, out)[0]
        g = gnina_rescore(receptor, out, resc)
    except Exception:
        continue
    mm = Chem.RemoveHs(m)
    rows.append({"chembl_id": r.id, "smiles": r.smiles, "active_prob": getattr(r, "active_prob", np.nan),
                 "MW": Descriptors.MolWt(mm), "LogP": Crippen.MolLogP(mm),
                 "TPSA": rdMolDescriptors.CalcTPSA(mm), "HBD": rdMolDescriptors.CalcNumHBD(mm),
                 "HBA": rdMolDescriptors.CalcNumHBA(mm), "RotB": rdMolDescriptors.CalcNumRotatableBonds(mm),
                 "smina_affinity": float(props(best)["minimizedAffinity"]),
                 "gnina_cnn_score": float(props(g)["CNNscore"]),
                 "gnina_cnn_affinity": float(props(g)["CNNaffinity"])})

res = pd.DataFrame(rows).sort_values("smina_affinity").reset_index(drop=True)
res["rank_smina"] = res["smina_affinity"].rank(method="min")
res["rank_cnn"] = res["gnina_cnn_score"].rank(ascending=False, method="min")
top10 = res.head(10).copy()
top10["combined_rank"] = top10["rank_smina"] + top10["rank_cnn"]
hits = top10.sort_values(["combined_rank", "smina_affinity"]).head(3)
res.to_csv(DOCK / "docking_results_all.csv", index=False)
hits.to_csv(DOCK / "table5_top3_hits.csv", index=False)
print(hits[["chembl_id", "active_prob", "smina_affinity", "gnina_cnn_affinity"]])
