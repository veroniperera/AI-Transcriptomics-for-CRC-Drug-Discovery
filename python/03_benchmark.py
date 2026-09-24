import json

import numpy as np
import optuna
import pandas as pd
from sklearn.model_selection import StratifiedKFold

import models
from common import ML, PROC, SEED, seed_all
from metrics import early_enrichment, metric_row
from features import DESC_2D, DESC_3D

N_TRIALS = 6
FEATS = DESC_2D + DESC_3D
seed_all()
optuna.logging.set_verbosity(optuna.logging.WARNING)

df = pd.read_csv(PROC / "features.csv")
tr, te = df[df.split == "train"], df[df.split == "test"]
Xtr, ytr = tr[FEATS].values, tr["label"].values
Xte, yte = te[FEATS].values, te["label"].values
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
folds = list(skf.split(Xtr, ytr))


def cv_metrics(name, params):
    rows = []
    for a, b in folds:
        m = models.build(name, params).fit(Xtr[a], ytr[a])
        rows.append(metric_row(ytr[b], m.predict_proba(Xtr[b])[:, 1]))
    return pd.DataFrame(rows)


cv_rows, hold_rows, best_params, preds = [], [], {}, {}
for name in models.NAMES:
    def objective(trial):
        params = models.suggest(name, trial)
        cv = cv_metrics(name, params)
        trial.set_user_attr("cv_mean", cv.mean().to_dict())
        trial.set_user_attr("cv_std", cv.std().to_dict())
        return cv["PR_AUC"].mean()

    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=SEED))
    study.optimize(objective, n_trials=N_TRIALS)
    bt = study.best_trial
    params = models.suggest(name, optuna.trial.FixedTrial(bt.params))
    best_params[name] = params

    cv_rows.append({"model": name, **{f"{k}_mean": v for k, v in bt.user_attrs["cv_mean"].items()},
                    **{f"{k}_sd": v for k, v in bt.user_attrs["cv_std"].items()}})

    fit = models.build(name, params).fit(Xtr, ytr)
    p = fit.predict_proba(Xte)[:, 1]
    preds[name] = p
    hold_rows.append({"model": name, **metric_row(yte, p)})
    print(name, round(bt.value, 3), round(hold_rows[-1]["PR_AUC"], 3))

cv_tab = pd.DataFrame(cv_rows).sort_values("PR_AUC_mean", ascending=False)
hold_tab = pd.DataFrame(hold_rows).sort_values("PR_AUC", ascending=False)
cv_tab.to_csv(ML / "table2_cv_metrics.csv", index=False)
hold_tab.to_csv(ML / "table3_holdout_metrics.csv", index=False)
json.dump(best_params, open(ML / "best_params.json", "w"), indent=2)

pd.DataFrame({"chembl_id": te["chembl_id"].values, "label": yte, **preds}).to_csv(
    ML / "holdout_predictions.csv", index=False)

best = hold_tab.iloc[0]["model"]
pd.DataFrame(early_enrichment(yte, preds[best])).to_csv(ML / "fig9_early_enrichment.csv", index=False)
open(ML / "best_model.txt", "w").write(best)
print("best:", best)
