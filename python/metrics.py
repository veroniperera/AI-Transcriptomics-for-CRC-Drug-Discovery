import numpy as np
from sklearn.metrics import average_precision_score, f1_score, matthews_corrcoef, roc_auc_score


def metric_row(y, p, thr=0.5):
    y = np.asarray(y)
    p = np.asarray(p)
    pred = (p >= thr).astype(int)
    return {"PR_AUC": average_precision_score(y, p), "ROC_AUC": roc_auc_score(y, p),
            "MCC": matthews_corrcoef(y, pred), "F1": f1_score(y, pred)}


def early_enrichment(y, p, fractions=(0.01, 0.02, 0.05, 0.10, 0.15, 0.20)):
    y = np.asarray(y)
    order = np.argsort(-np.asarray(p), kind="stable")
    base = y.sum() / len(y)
    rows = []
    for f in fractions:
        k = max(1, int(np.floor(f * len(y))))
        hits = int(y[order[:k]].sum())
        rows.append({"fraction": f, "k": k, "hits": hits, "n_active": int(y.sum()),
                     "n_total": len(y), "EF": (hits / k) / base})
    return rows
