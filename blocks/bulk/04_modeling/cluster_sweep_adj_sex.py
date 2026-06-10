#!/usr/bin/env python
"""
Cluster sweep adjusted for sex (and age): for each cluster Cx, fit
    qst_ppt_tr_avg_v1 ~ Cx + age + sex_male
and report Cx's adjusted coefficient, p-value, 95% CI, and partial correlation.
Question: does ANY single cluster predict pain after accounting for sex?

Also reports the raw (unadjusted) Cx-vs-pain correlation for contrast, and
multiple-testing context (6 clusters tested -> Bonferroni alpha = 0.05/6).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf

from regression_qst_v1 import load_data

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "regression_qst_v1_C2_out")  # keep sweep with the focused-model outputs
os.makedirs(OUTDIR, exist_ok=True)

TARGET = "qst_ppt_tr_avg_v1"
CLUSTERS = ["C1", "C2", "C3", "C4", "C5", "C6"]
N_TESTS = len(CLUSTERS)
BONF = 0.05 / N_TESTS


def main():
    X_all, y = load_data()
    df = X_all[CLUSTERS + ["age", "sex_male"]].copy()
    df[TARGET] = y.values

    rows = []
    for c in CLUSTERS:
        # raw (unadjusted) correlation
        sub = df[[c, TARGET]].dropna()
        r_raw, p_raw = stats.pearsonr(sub[c], sub[TARGET])

        # adjusted model: Cx + age + sex
        m = smf.ols(f"Q('{TARGET}') ~ {c} + age + sex_male", data=df, missing="drop").fit()
        ci = m.conf_int().loc[c]

        # partial correlation of Cx with pain, removing age + sex
        ok = df[[c, "age", "sex_male", TARGET]].dropna()
        ry = smf.ols(f"Q('{TARGET}') ~ age + sex_male", data=ok).fit().resid
        rc = smf.ols(f"{c} ~ age + sex_male", data=ok).fit().resid
        r_par, p_par = stats.pearsonr(rc, ry)

        rows.append({
            "cluster": c,
            "raw_r": r_raw, "raw_p": p_raw,
            "adj_coef": m.params[c], "adj_p": m.pvalues[c],
            "adj_CI_low": ci[0], "adj_CI_high": ci[1],
            "partial_r": r_par,
            "sig_raw": p_raw < 0.05,
            "sig_adj": m.pvalues[c] < 0.05,
            "sig_adj_bonf": m.pvalues[c] < BONF,
        })

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUTDIR, "cluster_sweep_adj_sex.csv"), index=False)

    show = res[["cluster", "raw_r", "raw_p", "adj_coef", "adj_p", "partial_r"]].copy()
    print(f"n = (per model, after dropping missing age) ; tested {N_TESTS} clusters; "
          f"Bonferroni alpha = {BONF:.4f}\n")
    print(show.round(4).to_string(index=False))
    n_raw = int(res["sig_raw"].sum()); n_adj = int(res["sig_adj"].sum())
    print(f"\nSignificant (p<0.05)  RAW: {n_raw}/{N_TESTS}  ->  ADJUSTED for age+sex: {n_adj}/{N_TESTS}")
    print(f"Surviving Bonferroni (p<{BONF:.4f}) after adjustment: {int(res['sig_adj_bonf'].sum())}/{N_TESTS}")

    # -------- figure: adjusted coefficients with 95% CI, colored by significance -------- #
    fig, ax = plt.subplots(figsize=(8.5, 5))
    order = res.iloc[::-1]  # C1 at bottom
    ypos = np.arange(len(order))
    err = np.vstack([order["adj_coef"] - order["adj_CI_low"],
                     order["adj_CI_high"] - order["adj_coef"]])
    colors = ["#55A868" if s else "#9aa0a6" for s in order["sig_adj"]]
    ax.errorbar(order["adj_coef"], ypos, xerr=err, fmt="o", ms=8,
                ecolor="gray", elinewidth=1.5, capsize=4,
                mfc="none", mec="none", zorder=1)
    ax.scatter(order["adj_coef"], ypos, s=80, c=colors, zorder=2, edgecolor="k", linewidth=.4)
    ax.axvline(0, color="r", ls="--", lw=1)
    ax.set_yticks(ypos); ax.set_yticklabels(order["cluster"])
    ax.set_xlabel("adjusted coefficient for cluster  (qst ~ Cx + age + sex)")
    ax.set_title("Each cluster's effect on pain, adjusted for age + sex\n"
                 "(green = p<0.05; gray = not significant; line = 95% CI)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "cluster_sweep_adj_sex.png"), dpi=160)
    plt.close()

    print(f"\nOutputs -> {OUTDIR}/cluster_sweep_adj_sex.csv , cluster_sweep_adj_sex.png")


if __name__ == "__main__":
    main()
