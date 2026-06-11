#!/usr/bin/env python
"""
Compare regression models predicting visit-1 pressure-pain threshold
(qst_ppt_tr_avg_v1) from per-sample k=6 cluster mean z-scores (age/sex excluded).

Design notes (n ~ 49, p = 6 -> small-n regime):
  - Lean, regularized model set; ensembles deliberately excluded (overfit at n=49).
  - Nested CV: hyperparameters tuned in an inner loop so reported scores are honest.
  - All preprocessing (impute age, standardize) lives INSIDE the pipeline, so it is
    fit on the training fold only -> no leakage from the held-out fold.
  - Pure regression: ROC-AUC is undefined for a continuous target, so we report
    RMSE/MAE/R2/Spearman/Pearson and an out-of-fold predicted-vs-actual plot instead.

Inputs:
  kmeans_k6_cluster_mean_zscore_per_sample_LFC_P5P95_ZSCORE.csv  (long: Cluster,SampleID,SampleLab,MeanZScore)
  example_pain_metadata.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.base import clone
from sklearn.inspection import permutation_importance
from sklearn.model_selection import (
    RepeatedKFold, KFold, LeaveOneOut, GridSearchCV,
    cross_validate, cross_val_predict,
)
from sklearn.metrics import mean_squared_error, mean_absolute_error, make_scorer

warnings.filterwarnings("ignore")  # silence convergence chatter on tiny folds
RNG = 42

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV = os.path.join(HERE, "kmeans_k6_cluster_mean_zscore_per_sample_LFC_P5P95_ZSCORE.csv")
META_CSV  = os.path.join(HERE, "example_pain_metadata.csv")
OUTDIR    = os.path.join(HERE, "phenotype_regression_nestedCV_out")
os.makedirs(OUTDIR, exist_ok=True)

TARGET   = "qst_ppt_tr_avg_v1"
CLUSTERS = ["C1", "C2", "C3", "C4", "C5", "C6"]
FEATURES = CLUSTERS  # age/sex deliberately excluded: clusters-only predictors


# --------------------------------------------------------------------------- #
# 1. Load + assemble the modeling table
# --------------------------------------------------------------------------- #
def normalize_id(sample_lab: str) -> str:
    """'PROMO14_V1' -> 'PROMO_14' to match metadata study_id (keeps zero-padding)."""
    import re
    m = re.match(r"PROMO(\d+)_V\d+", str(sample_lab))
    return f"PROMO_{m.group(1)}" if m else str(sample_lab)


def load_data():
    long = pd.read_csv(INPUT_CSV)
    wide = (long.pivot(index="SampleLab", columns="Cluster", values="MeanZScore")
                .reset_index())

    # Save a wide-format copy of the input CSV (one row per sample, one column per cluster).
    wide_csv = os.path.join(HERE, "kmeans_k6_cluster_mean_zscore_per_sample_LFC_P5P95_ZSCORE_wide.csv")
    wide.to_csv(wide_csv, index=False)
    print(f"[wide ] wrote wide-format copy -> {wide_csv}")

    wide["study_id"] = wide["SampleLab"].map(normalize_id)

    # Collapse technical replicates (e.g. PROMO24_V1 / _V1.1 / _V1.2 -> one patient):
    # average each cluster's z-score so every patient contributes exactly one row.
    rep = wide.groupby("study_id").size()
    if (rep > 1).any():
        print(f"[reps ] averaging replicates for: {rep[rep > 1].to_dict()}")
    wide = wide.groupby("study_id", as_index=False)[CLUSTERS].mean()

    # age/sex are intentionally excluded from the predictors; only pull the target.
    meta = pd.read_csv(META_CSV)
    keep = ["study_id", TARGET]
    meta = meta[keep]

    df = wide.merge(meta, on="study_id", how="inner")

    print(f"[merge] input samples={wide.shape[0]}, metadata rows={meta.shape[0]}, "
          f"merged={df.shape[0]}")
    unmatched = set(wide['study_id']) ^ set(meta['study_id'])
    if unmatched:
        print(f"[merge] WARNING unmatched ids ({len(unmatched)}): {sorted(unmatched)}")

    # Require a non-missing target; cluster predictors are complete per sample.
    df = df.dropna(subset=[TARGET]).reset_index(drop=True)

    print(f"[data ] n={df.shape[0]} with non-missing target; "
          f"predictors = {len(FEATURES)} clusters (age/sex excluded)")

    X = df[FEATURES].copy()
    y = df[TARGET].copy()
    ids = df["study_id"].copy()
    return X, y, ids


# --------------------------------------------------------------------------- #
# 2. Collinearity check (VIF) on the predictors
# --------------------------------------------------------------------------- #
def vif_table(X: pd.DataFrame) -> pd.DataFrame:
    Xz = (X - X.mean()) / X.std(ddof=0)
    Xz = Xz.fillna(0.0)  # the single missing age -> mean (0) just for the diagnostic
    vifs = []
    for j, col in enumerate(Xz.columns):
        others = [c for c in Xz.columns if c != col]
        r2 = LinearRegression().fit(Xz[others], Xz[col]).score(Xz[others], Xz[col])
        vifs.append(1.0 / (1.0 - r2) if r2 < 1 else np.inf)
    out = pd.DataFrame({"feature": Xz.columns, "VIF": np.round(vifs, 3)})
    out.to_csv(os.path.join(OUTDIR, "vif.csv"), index=False)
    print("\n[VIF] (>5 = notable collinearity, >10 = severe)")
    print(out.to_string(index=False))
    return out


# --------------------------------------------------------------------------- #
# 3. Model zoo: each = (pipeline, param_grid for inner tuning)
# --------------------------------------------------------------------------- #
def make_pipe(estimator):
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", estimator),
    ])


def model_zoo():
    return {
        "OLS":             (make_pipe(LinearRegression()), {}),
        "Ridge":           (make_pipe(Ridge(random_state=RNG)),
                            {"model__alpha": np.logspace(-2, 3, 12)}),
        "Lasso":           (make_pipe(Lasso(max_iter=20000, random_state=RNG)),
                            {"model__alpha": np.logspace(-3, 1, 12)}),
        "ElasticNet":      (make_pipe(ElasticNet(max_iter=20000, random_state=RNG)),
                            {"model__alpha": np.logspace(-3, 1, 8),
                             "model__l1_ratio": [0.2, 0.5, 0.8]}),
        "SVR (RBF)":       (make_pipe(SVR(kernel="rbf")),
                            {"model__C": [0.1, 1, 10, 100],
                             "model__gamma": ["scale", 0.01, 0.1]}),
        "kNN":             (make_pipe(KNeighborsRegressor()),
                            {"model__n_neighbors": [3, 5, 7, 9, 11]}),
    }


# --------------------------------------------------------------------------- #
# 4. Nested CV scoring + out-of-fold predictions
# --------------------------------------------------------------------------- #
def rmse(y, yhat):
    return np.sqrt(mean_squared_error(y, yhat))


def _pearson_r2(y_true, y_pred):
    """Squared Pearson r — always in [0, 1]. Returns 0 when predictions are constant."""
    r, _ = stats.pearsonr(y_true, y_pred)
    return 0.0 if np.isnan(r) else r ** 2


def evaluate(X, y):
    outer = RepeatedKFold(n_splits=5, n_repeats=10, random_state=RNG)
    inner = KFold(n_splits=5, shuffle=True, random_state=RNG)
    loo = LeaveOneOut()
    scoring = {
        "rmse": "neg_root_mean_squared_error",
        "mae": "neg_mean_absolute_error",
        # Squared Pearson r computed per fold, then averaged — stays in [0, 1].
        "r2": make_scorer(_pearson_r2),
    }

    rows, oof_preds = [], {}
    for name, (pipe, grid) in model_zoo().items():
        est = (GridSearchCV(pipe, grid, cv=inner,
                            scoring="neg_root_mean_squared_error", n_jobs=-1)
               if grid else pipe)

        cv = cross_validate(est, X, y, cv=outer, scoring=scoring, n_jobs=-1)
        # LOO out-of-fold predictions (one prediction per sample) for scatter plots only.
        yhat = cross_val_predict(est, X, y, cv=loo, n_jobs=-1)
        oof_preds[name] = yhat

        rho = stats.spearmanr(y, yhat).correlation
        r_p = stats.pearsonr(y, yhat)[0]
        rows.append({
            "model": name,
            # In-fold squared-Pearson R² averaged across all 5x10 nested-CV folds.
            "R2": cv["test_r2"].mean(),
            "R2_std": cv["test_r2"].std(),
            "RMSE_mean": -cv["test_rmse"].mean(),  "RMSE_std": cv["test_rmse"].std(),
            "MAE_mean": -cv["test_mae"].mean(),    "MAE_std": cv["test_mae"].std(),
            "Spearman_oof": rho, "Pearson_oof": r_p,
            "RMSE_oof": rmse(y, yhat),
        })
        print(f"[cv] {name:16s}  R2(in-fold)={rows[-1]['R2']:.3f}  "
              f"RMSE={rows[-1]['RMSE_mean']:.3f}  Spearman={rho:+.3f}")

    res = pd.DataFrame(rows).sort_values("RMSE_mean").reset_index(drop=True)
    res.to_csv(os.path.join(OUTDIR, "results_table.csv"), index=False)
    return res, oof_preds


# --------------------------------------------------------------------------- #
# 5. Plots
# --------------------------------------------------------------------------- #
def plot_comparison(res):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    r = res.sort_values("RMSE_mean")
    axes[0].barh(r["model"], r["RMSE_mean"], xerr=r["RMSE_std"],
                 color="#4C72B0", alpha=.85)
    axes[0].set_xlabel("CV RMSE (lower = better)"); axes[0].invert_yaxis()
    axes[0].set_title("Cross-validated RMSE (5x10 nested CV)")

    r2 = res.sort_values("R2", ascending=True)
    axes[1].barh(r2["model"], r2["R2"], color="#55A868", alpha=.85)
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("R$^2$ (squared Pearson r, in-fold mean, 0–1)")
    axes[1].set_title("In-fold CV R$^2$ (5×10 repeated nested KFold)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "model_comparison.png"), dpi=160)
    plt.close()


def plot_oof(res, oof_preds, y):
    # study_id labels (same row order as X/y -> align by position) for outlier annotation.
    ids = pd.read_csv(os.path.join(OUTDIR, "modeling_table.csv"))["study_id"].tolist()
    yv = y.values
    models = res["model"].tolist()
    r2_lookup = res.set_index("model")["R2"]
    n = len(models)
    ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow), squeeze=False)
    lo, hi = y.min(), y.max()
    for i, name in enumerate(models):
        ax = axes[i // ncol][i % ncol]
        yhat = oof_preds[name]
        ax.scatter(yv, yhat, s=28, alpha=.7, edgecolor="k", linewidth=.3)
        ax.plot([lo, hi], [lo, hi], "r--", lw=1)
        # Flag + label outliers: residual (distance from the y=x line) > 2 SD for this model.
        resid = yv - yhat
        thr = 2 * resid.std()
        out = np.where(np.abs(resid) > thr)[0]
        ax.scatter(yv[out], yhat[out], s=46, facecolors="none",
                   edgecolors="#C44E52", linewidth=1.4, zorder=3)
        for k in out:
            ax.annotate(ids[k], (yv[k], yhat[k]), fontsize=7, color="#C44E52",
                        textcoords="offset points", xytext=(4, 3))
        r2 = r2_lookup[name]  # in-fold CV R² from nested cross-validation
        ax.set_title(f"{name}\nCV R$^2$={r2:.2f}", fontsize=10)
        ax.set_xlabel("actual (LOO)"); ax.set_ylabel("predicted (LOO)")
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("LOO predicted vs actual — qst_ppt_tr_avg_v1 (PROMO12 excluded)", y=1.0)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "pred_vs_actual.png"), dpi=160)
    plt.close()


def plot_residuals(res, oof_preds, y):
    best = res.iloc[0]["model"]
    yhat = oof_preds[best]
    resid = y.values - yhat
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].scatter(yhat, resid, s=28, alpha=.7, edgecolor="k", linewidth=.3)
    axes[0].axhline(0, color="r", ls="--", lw=1)
    axes[0].set_xlabel("predicted"); axes[0].set_ylabel("residual")
    axes[0].set_title(f"Residuals vs fitted — {best}")
    stats.probplot(resid, plot=axes[1])
    axes[1].set_title("Residual Q-Q")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "residuals.png"), dpi=160)
    plt.close()


def plot_coefficients(X, y):
    """Standardized coefficients from the linear models (fit on full standardized data)."""
    fitted = {}
    for name, est in [("OLS", LinearRegression()),
                      ("Ridge", Ridge(alpha=1.0)),
                      ("Lasso", Lasso(alpha=0.05, max_iter=20000)),
                      ("ElasticNet", ElasticNet(alpha=0.05, l1_ratio=0.5, max_iter=20000))]:
        pipe = make_pipe(est).fit(X, y)
        fitted[name] = pipe.named_steps["model"].coef_
    coef = pd.DataFrame(fitted, index=FEATURES)
    coef.to_csv(os.path.join(OUTDIR, "linear_coefficients.csv"))

    ax = coef.plot(kind="bar", figsize=(11, 5), width=.8)
    ax.axhline(0, color="k", lw=.8)
    ax.set_ylabel("standardized coefficient")
    ax.set_title("Linear-model coefficients (sign/size of each predictor's effect)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "coefficients.png"), dpi=160)
    plt.close()


def perm_importance(X, y):
    """Held-out permutation importance: for each model, refit on each training
    fold, then permute one cluster at a time on the TEST fold and measure how
    much RMSE rises. Positive = the cluster helps prediction; ~0 = it doesn't.
    Averaged across 5 folds. (scoring = RMSE, not the default R2, which is too
    noisy on ~9-sample test folds.)"""
    cv = KFold(n_splits=5, shuffle=True, random_state=RNG)
    inner = KFold(n_splits=5, shuffle=True, random_state=RNG)
    cols = {}
    for name, (pipe, grid) in model_zoo().items():
        fold_imps = []
        for tr, te in cv.split(X):
            est = (GridSearchCV(clone(pipe), grid, cv=inner,
                                scoring="neg_root_mean_squared_error", n_jobs=-1)
                   if grid else clone(pipe))
            est.fit(X.iloc[tr], y.iloc[tr])
            r = permutation_importance(est, X.iloc[te], y.iloc[te], n_repeats=5,
                                       random_state=RNG,
                                       scoring="neg_root_mean_squared_error")
            fold_imps.append(r.importances_mean)
        cols[name] = np.mean(fold_imps, axis=0)

    imp = pd.DataFrame(cols, index=FEATURES)
    imp.to_csv(os.path.join(OUTDIR, "permutation_importance.csv"))

    ax = imp.plot(kind="bar", figsize=(11, 5), width=.8)
    ax.axhline(0, color="k", lw=.8)
    ax.set_ylabel("RMSE increase when cluster is permuted")
    ax.set_title("Permutation importance (held-out): how much each cluster helps prediction "
                 "(higher = more useful; ~0 = no contribution)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "permutation_importance.png"), dpi=160)
    plt.close()

    print("\n[perm] held-out permutation importance (RMSE rise when permuted; ~0 = no signal):")
    print(imp.round(4).to_string())
    return imp


def plot_promo12_comparison(r2_compare):
    models = r2_compare["model"].tolist()
    x = np.arange(len(models))
    w = 0.35

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: grouped bars — with vs without PROMO12
    ax = axes[0]
    bars_with    = ax.bar(x - w/2, r2_compare["R2_with_PROMO12"],    w,
                          label="With PROMO12",    color="#4C72B0", alpha=.85)
    bars_without = ax.bar(x + w/2, r2_compare["R2_without_PROMO12"], w,
                          label="Without PROMO12", color="#55A868", alpha=.85)
    ax.set_xticks(x); ax.set_xticklabels(models, rotation=20, ha="right")
    ax.set_ylim(0, min(1.0, r2_compare["R2_without_PROMO12"].max() * 1.35))
    ax.set_ylabel("R$^2$ (in-fold squared Pearson, 0–1)")
    ax.set_title("In-fold CV R$^2$: with vs without PROMO12")
    ax.legend()
    for bar in list(bars_with) + list(bars_without):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.004,
                f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=8)

    # Right: delta (R² change from removing PROMO12)
    ax2 = axes[1]
    colors = ["#C44E52" if d < 0 else "#55A868" for d in r2_compare["R2_change"]]
    ax2.bar(x, r2_compare["R2_change"], color=colors, alpha=.85)
    ax2.axhline(0, color="k", lw=0.8)
    ax2.set_xticks(x); ax2.set_xticklabels(models, rotation=20, ha="right")
    ax2.set_ylabel("ΔR$^2$ (without − with PROMO12)")
    ax2.set_title("R$^2$ change after removing PROMO12\n(green = improves, red = worsens)")
    for i, (xi, d) in enumerate(zip(x, r2_compare["R2_change"])):
        ax2.text(xi, d + (0.002 if d >= 0 else -0.004),
                 f"{d:+.3f}", ha="center", va="bottom" if d >= 0 else "top", fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "promo12_r2_comparison.png"), dpi=160)
    plt.close()


# --------------------------------------------------------------------------- #
def main():
    X_full, y_full, ids_full = load_data()

    # --- R² with PROMO12 included ---
    print("\n--- Evaluating WITH PROMO12 (n={}) ---".format(len(y_full)))
    res_with, _ = evaluate(X_full, y_full)

    # --- Remove PROMO12 ---
    mask = ids_full != "PROMO_12"
    X = X_full[mask].reset_index(drop=True)
    y = y_full[mask].reset_index(drop=True)
    ids = ids_full[mask].reset_index(drop=True)
    print(f"\n--- Evaluating WITHOUT PROMO12 (n={len(y)}) ---")
    res, oof = evaluate(X, y)

    # Write modeling table for the primary (PROMO12-excluded) analysis
    pd.concat([ids, X, y], axis=1).to_csv(
        os.path.join(OUTDIR, "modeling_table.csv"), index=False)

    # --- Side-by-side R² comparison ---
    r2_compare = (res_with[["model", "R2"]].rename(columns={"R2": "R2_with_PROMO12"})
                  .merge(res[["model", "R2"]].rename(columns={"R2": "R2_without_PROMO12"}),
                         on="model"))
    r2_compare["R2_change"] = r2_compare["R2_without_PROMO12"] - r2_compare["R2_with_PROMO12"]
    r2_compare.to_csv(os.path.join(OUTDIR, "r2_promo12_comparison.csv"), index=False)
    plot_promo12_comparison(r2_compare)
    print("\n===== R² COMPARISON (in-fold squared Pearson) =====")
    print(r2_compare.round(3).to_string(index=False))

    print("\n===== RESULTS WITHOUT PROMO12 (sorted by CV RMSE) =====")
    show = res[["model", "R2", "RMSE_mean", "MAE_mean", "Spearman_oof", "Pearson_oof"]]
    print(show.round(3).to_string(index=False))

    vif_table(X)
    plot_comparison(res)
    plot_oof(res, oof, y)
    plot_residuals(res, oof, y)
    plot_coefficients(X, y)
    perm_importance(X, y)
    print(f"\nAll outputs written to: {OUTDIR}")


if __name__ == "__main__":
    main()
