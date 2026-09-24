import json

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

import models
from common import ML, PROC, SEED, seed_all
from features import DESC_2D, DESC_3D, ad_score
from metrics import metric_row

N_PERM = 30
FEATS = DESC_2D + DESC_3D
seed_all()

df = pd.read_csv(PROC / "features.csv")
fp = np.load(PROC / "fingerprints.npz")["ecfp4"]
tr_mask, te_mask = (df.split == "train").values, (df.split == "test").values
Xtr, ytr = df.loc[tr_mask, FEATS].values, df.loc[tr_mask, "label"].values
Xte, yte = df.loc[te_mask, FEATS].values, df.loc[te_mask, "label"].values

name = open(ML / "best_model.txt").read().strip()
params = json.load(open(ML / "best_params.json"))[name]


def oof_metrics(y):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = np.zeros(len(y))
    for a, b in skf.split(Xtr, y):
        oof[b] = models.build(name, params).fit(Xtr[a], y[a]).predict_proba(Xtr[b])[:, 1]
    return metric_row(y, oof)


real = oof_metrics(ytr)
rng = np.random.RandomState(SEED)
null = pd.DataFrame([oof_metrics(rng.permutation(ytr)) for _ in range(N_PERM)])

rows = []
for m in real:
    mu, sd = null[m].mean(), null[m].std()
    p = (1 + (null[m] >= real[m]).sum()) / (1 + N_PERM)
    rows.append({"metric": m, "real": real[m], "rand_mean": mu, "rand_sd": sd,
                 "z": (real[m] - mu) / sd, "p": p, "significant": p < 0.05})
pd.DataFrame(rows).to_csv(ML / "table4_y_randomization.csv", index=False)
null.to_csv(ML / "y_randomization_null.csv", index=False)

final = models.build(name, params).fit(Xtr, ytr)
joblib.dump(final, ML / "final_model.joblib")
pte = final.predict_proba(Xte)[:, 1]

train_ad = ad_score(fp[tr_mask], fp[tr_mask], k=5, exclude_self=True)
test_ad = ad_score(fp[te_mask], fp[tr_mask], k=5)
thr_p5 = float(np.percentile(train_ad, 5))
thr_med = float(np.median(test_ad))
json.dump({"train_mean": float(train_ad.mean()), "test_mean": float(test_ad.mean()),
           "threshold_p5_train": thr_p5, "n_out_domain_at_p5": int((test_ad < thr_p5).sum()),
           "threshold_median_test": thr_med, "n_in_domain": int((test_ad >= thr_med).sum()),
           "n_out_domain": int((test_ad < thr_med).sum())},
          open(ML / "applicability_domain.json", "w"), indent=2)

rows = []
for label, mask in (("in_domain", test_ad >= thr_med), ("out_domain", test_ad < thr_med)):
    if len(np.unique(yte[mask])) == 2:
        rows.append({"domain": label, "n": int(mask.sum()), **metric_row(yte[mask], pte[mask])})
pd.DataFrame(rows).to_csv(ML / "ad_domain_metrics.csv", index=False)

cov = []
for q in np.linspace(0, 95, 20):
    t = np.percentile(test_ad, q)
    mask = test_ad >= t
    if mask.sum() > 4 and len(np.unique(yte[mask])) == 2:
        cov.append({"percentile": q, "threshold": t, "coverage": mask.mean(),
                    "PR_AUC": metric_row(yte[mask], pte[mask])["PR_AUC"]})
pd.DataFrame(cov).to_csv(ML / "ad_coverage_curve.csv", index=False)
print(pd.DataFrame(rows).round(3))
