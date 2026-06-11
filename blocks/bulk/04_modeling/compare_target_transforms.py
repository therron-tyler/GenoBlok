#!/usr/bin/env python
"""Compare target transforms for the clusters-only PPT regression.

For each model in the zoo, evaluate it under three target settings:
  - none         : predict PPT directly (current setup)
  - log          : predict log(PPT), invert back to PPT for scoring
  - yeo-johnson  : PowerTransformer fit per-fold, invert back for scoring

All transforms are applied INSIDE the CV via TransformedTargetRegressor, so
the transform (and Yeo-Johnson's lambda) is fit on training folds only -> no
leakage. Scores are computed on the ORIGINAL PPT scale, so they are directly
comparable to results_table.csv.

R2 here = squared Pearson correlation of LOO out-of-fold preds vs actual (0-1).
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

from sklearn.compose import TransformedTargetRegressor
from sklearn.preprocessing import PowerTransformer
from sklearn.model_selection import (
    RepeatedKFold, KFold, LeaveOneOut, GridSearchCV,
    cross_validate, cross_val_predict,
)
from sklearn.metrics import mean_squared_error

from phenotype_regression_nestedCV import load_data, model_zoo

RNG = 42
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "phenotype_regression_nestedCV_out")


def wrap_target(est, transform):
    if transform == "none":
        return est
    if transform == "log":
        return TransformedTargetRegressor(regressor=est, func=np.log, inverse_func=np.exp)
    if transform == "yeo-johnson":
        return TransformedTargetRegressor(
            regressor=est, transformer=PowerTransformer(method="yeo-johnson"))
    raise ValueError(transform)


def main():
    X, y = load_data()

    outer = RepeatedKFold(n_splits=5, n_repeats=10, random_state=RNG)
    inner = KFold(n_splits=5, shuffle=True, random_state=RNG)
    loo = LeaveOneOut()

    rows = []
    for transform in ["none", "log", "yeo-johnson"]:
        for name, (pipe, grid) in model_zoo().items():
            base = (GridSearchCV(pipe, grid, cv=inner,
                                 scoring="neg_root_mean_squared_error", n_jobs=-1)
                    if grid else pipe)
            est = wrap_target(base, transform)

            cv = cross_validate(est, X, y, cv=outer,
                                scoring="neg_root_mean_squared_error", n_jobs=-1)
            yhat = cross_val_predict(est, X, y, cv=loo, n_jobs=-1)

            r_p = stats.pearsonr(y, yhat)[0]
            rho = stats.spearmanr(y, yhat).correlation
            rows.append({
                "transform": transform,
                "model": name,
                "R2": r_p ** 2,                 # raw squared-Pearson, original PPT scale
                "RMSE": -cv["test_score"].mean(),
                "Pearson": r_p,
                "Spearman": rho,
            })
            print(f"[{transform:11s}] {name:12s}  R2={r_p**2:.3f}  "
                  f"RMSE={-cv['test_score'].mean():.3f}  Pearson={r_p:+.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUTDIR, "target_transform_comparison.csv"), index=False)

    # Tidy side-by-side pivots for quick reading.
    print("\n===== R2 (squared Pearson, original PPT scale) =====")
    r2p = res.pivot(index="model", columns="transform", values="R2")[["none", "log", "yeo-johnson"]]
    print(r2p.round(3).to_string())

    print("\n===== RMSE (original PPT scale, lower=better) =====")
    rmp = res.pivot(index="model", columns="transform", values="RMSE")[["none", "log", "yeo-johnson"]]
    print(rmp.round(3).to_string())

    # Best cell per metric, for a one-line takeaway.
    best_r2 = res.loc[res["R2"].idxmax()]
    best_rmse = res.loc[res["RMSE"].idxmin()]
    print(f"\nBest R2  : {best_r2['model']} + {best_r2['transform']}  -> R2={best_r2['R2']:.3f}")
    print(f"Best RMSE: {best_rmse['model']} + {best_rmse['transform']}  -> RMSE={best_rmse['RMSE']:.3f}")
    print(f"\nSaved: {os.path.join(OUTDIR, 'target_transform_comparison.csv')}")


if __name__ == "__main__":
    main()
