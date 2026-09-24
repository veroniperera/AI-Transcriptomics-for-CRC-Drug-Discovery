import argparse
import subprocess

import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.metrics import roc_auc_score

from common import BOX_CENTER, BOX_SIZE, DOCK, DUDE, SEED, seed_all
from features import embed

ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
d = sub.add_parser("dock")
d.add_argument("--actives", required=True)
d.add_argument("--decoys", required=True)
d.add_argument("--receptor", default=str(DOCK / "work" / "receptor.pdbqt"))
d.add_argument("--n_actives", type=int, default=10)
d.add_argument("--ratio", type=int, default=50)
d.add_argument("--cpu", type=int, default=8)
a = sub.add_parser("analyze")
a.add_argument("--csv", default=str(DUDE / "dude_docking_SMALLBOX.csv"))
a.add_argument("--n_boot", type=int, default=2000)
args = ap.parse_args()
seed_all()


def read_ism(path):
    rows = [l.split() for l in open(path) if l.strip()]
    return pd.DataFrame({"ID": [r[1] if len(r) > 1 else f"cmpd{i}" for i, r in enumerate(rows)],
                         "Smiles": [r[0] for r in rows]})


def vina_score(receptor, smiles, name, work):
    m = embed(Chem.MolFromSmiles(smiles))
    if m is None:
        return np.nan
    sdf, pdbqt, out = work / f"{name}.sdf", work / f"{name}.pdbqt", work / f"{name}_out.pdbqt"
    w = Chem.SDWriter(str(sdf))
    w.write(m)
    w.close()
    try:
        subprocess.run(["obabel", str(sdf), "-O", str(pdbqt)], check=True, capture_output=True)
        subprocess.run(["vina", "--receptor", receptor, "--ligand", str(pdbqt),
                        "--center_x", str(BOX_CENTER[0]), "--center_y", str(BOX_CENTER[1]), "--center_z", str(BOX_CENTER[2]),
                        "--size_x", str(BOX_SIZE[0]), "--size_y", str(BOX_SIZE[1]), "--size_z", str(BOX_SIZE[2]),
                        "--exhaustiveness", "8", "--cpu", str(args.cpu), "--seed", str(SEED),
                        "--out", str(out)], check=True, capture_output=True)
        for line in open(out):
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    except Exception:
        return np.nan
    return np.nan


def ef(y, s, frac):
    k = int(np.floor(frac * len(y)))
    top = np.argsort(-s, kind="stable")[:k]
    return (y[top].sum() / k) / (y.sum() / len(y))


if args.cmd == "dock":
    work = DUDE / "work"
    work.mkdir(exist_ok=True)
    act = read_ism(args.actives).sample(n=args.n_actives, random_state=SEED)
    dec = read_ism(args.decoys).sample(n=args.n_actives * args.ratio, random_state=SEED)
    act["Bioactivity_Class"], act["Label"] = "Active", 1
    dec["Bioactivity_Class"], dec["Label"] = "Decoy", 0
    df = pd.concat([dec, act], ignore_index=True).sample(frac=1, random_state=SEED).reset_index(drop=True)
    df["Vina_Affinity"] = [vina_score(args.receptor, r.Smiles, f"m{i}", work) for i, r in enumerate(df.itertuples())]
    df.to_csv(DUDE / "dude_docking_SMALLBOX.csv", index=False)
    print(len(df), int(df["Vina_Affinity"].isna().sum()))

if args.cmd == "analyze":
    df = pd.read_csv(args.csv)
    n_sub, a_sub = len(df), int(df["Label"].sum())
    df = df.dropna(subset=["Vina_Affinity"]).copy()
    df["HA"] = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in df["Smiles"]]
    df["LE"] = -df["Vina_Affinity"] / df["HA"]
    y = df["Label"].values
    fracs = [0.01, 0.05, 0.10, 0.20]
    rng = np.random.default_rng(SEED)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    out = []
    for score_name, s in (("LE", df["LE"].values), ("Vina_raw", -df["Vina_Affinity"].values)):
        boot = []
        for _ in range(args.n_boot):
            idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
            boot.append([roc_auc_score(y[idx], s[idx])] + [ef(y[idx], s[idx], f) for f in fracs])
        ci = np.percentile(np.array(boot), [2.5, 97.5], axis=0)
        out.append({"score": score_name, "metric": "ROC_AUC", "value": roc_auc_score(y, s),
                    "ci_low": ci[0, 0], "ci_high": ci[1, 0]})
        for i, f in enumerate(fracs):
            k = int(np.floor(f * len(y)))
            hits = int(y[np.argsort(-s, kind="stable")[:k]].sum())
            out.append({"score": score_name, "metric": f"EF{int(f*100)}%", "value": ef(y, s, f),
                        "ci_low": ci[0, i + 1], "ci_high": ci[1, i + 1], "hits": hits, "k": k})
    res = pd.DataFrame(out)
    res.to_csv(DUDE / "dude_enrichment_results.csv", index=False)
    print(f"submitted {n_sub}/{a_sub}; ranked {len(df)}/{int(y.sum())}")
    print(res.round(2).to_string(index=False))
    print(df.groupby("Label")[["Vina_Affinity", "HA", "LE"]].mean().round(3))
