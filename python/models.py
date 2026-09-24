from lightgbm import LGBMClassifier
from sklearn.ensemble import (ExtraTreesClassifier, GradientBoostingClassifier,
                              HistGradientBoostingClassifier, RandomForestClassifier)
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

from common import SEED

NAMES = ["SVM_RBF", "LightGBM", "KNN", "RandomForest", "GradBoost",
         "XGBoost", "HistGradBoost", "ExtraTrees", "LogisticReg", "MLP"]


def suggest(name, t):
    if name in ("RandomForest", "ExtraTrees"):
        return dict(n_estimators=t.suggest_int("n_estimators", 100, 600),
                    max_depth=t.suggest_int("max_depth", 3, 20),
                    min_samples_split=t.suggest_int("min_samples_split", 2, 20),
                    max_features=t.suggest_categorical("max_features", ["sqrt", "log2"]))
    if name in ("XGBoost", "LightGBM"):
        return dict(n_estimators=t.suggest_int("n_estimators", 100, 600),
                    learning_rate=t.suggest_float("learning_rate", 0.01, 0.3, log=True),
                    subsample=t.suggest_float("subsample", 0.5, 1.0),
                    colsample_bytree=t.suggest_float("colsample_bytree", 0.5, 1.0),
                    max_depth=t.suggest_int("max_depth", 3, 10),
                    reg_alpha=t.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
                    reg_lambda=t.suggest_float("reg_lambda", 1e-8, 10.0, log=True))
    if name == "HistGradBoost":
        return dict(max_iter=t.suggest_int("max_iter", 100, 600),
                    max_depth=t.suggest_int("max_depth", 3, 12),
                    learning_rate=t.suggest_float("learning_rate", 0.01, 0.3, log=True),
                    l2_regularization=t.suggest_float("l2_regularization", 1e-8, 10.0, log=True))
    if name == "GradBoost":
        return dict(n_estimators=t.suggest_int("n_estimators", 100, 600),
                    max_depth=t.suggest_int("max_depth", 2, 8),
                    learning_rate=t.suggest_float("learning_rate", 0.01, 0.3, log=True),
                    subsample=t.suggest_float("subsample", 0.5, 1.0),
                    min_samples_leaf=t.suggest_int("min_samples_leaf", 1, 10))
    if name == "SVM_RBF":
        return dict(C=t.suggest_float("C", 1e-2, 1e3, log=True),
                    gamma=t.suggest_float("gamma", 1e-4, 1e1, log=True))
    if name == "MLP":
        n_layers = t.suggest_int("n_layers", 1, 3)
        units = t.suggest_categorical("units", [32, 64, 128, 256])
        return dict(hidden_layer_sizes=[units] * n_layers,
                    learning_rate_init=t.suggest_float("learning_rate_init", 1e-4, 1e-2, log=True),
                    alpha=t.suggest_float("alpha", 1e-6, 1e-1, log=True))
    if name == "LogisticReg":
        return dict(C=t.suggest_float("C", 1e-3, 1e2, log=True),
                    penalty=t.suggest_categorical("penalty", ["l1", "l2"]))
    if name == "KNN":
        return dict(n_neighbors=t.suggest_int("n_neighbors", 3, 25),
                    weights=t.suggest_categorical("weights", ["uniform", "distance"]),
                    metric=t.suggest_categorical("metric", ["euclidean", "manhattan", "chebyshev"]))
    raise ValueError(name)


def build(name, p):
    p = dict(p)
    if name == "RandomForest":
        clf = RandomForestClassifier(random_state=SEED, n_jobs=-1, **p)
    elif name == "ExtraTrees":
        clf = ExtraTreesClassifier(random_state=SEED, n_jobs=-1, **p)
    elif name == "XGBoost":
        clf = XGBClassifier(random_state=SEED, n_jobs=-1, eval_metric="logloss", verbosity=0, **p)
    elif name == "LightGBM":
        clf = LGBMClassifier(random_state=SEED, n_jobs=-1, verbose=-1, **p)
    elif name == "HistGradBoost":
        clf = HistGradientBoostingClassifier(random_state=SEED, **p)
    elif name == "GradBoost":
        clf = GradientBoostingClassifier(random_state=SEED, **p)
    elif name == "SVM_RBF":
        clf = SVC(kernel="rbf", probability=True, random_state=SEED, **p)
    elif name == "MLP":
        clf = MLPClassifier(random_state=SEED, max_iter=1000, hidden_layer_sizes=tuple(p.pop("hidden_layer_sizes")), **p)
    elif name == "LogisticReg":
        clf = LogisticRegression(solver="liblinear", max_iter=5000, random_state=SEED, **p)
    elif name == "KNN":
        clf = KNeighborsClassifier(**p)
    else:
        raise ValueError(name)
    return Pipeline([("scaler", StandardScaler()), ("clf", clf)])
