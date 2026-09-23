#!/usr/bin/env python3
"""
correlate_perturbation_with_features.py

Answers the actual question fibril_perturbation_analysis.py's output was
built for: does curvature, or local_sustained_gap_nm, change with
proximity to each ultrastructural feature type?

For each feature type found in the CSV (dist_to_ER, dist_to_mito, etc.),
this does two complementary things:

1. SPEARMAN CORRELATION between distance-to-feature and each metric
   (curvature, local_sustained_gap_nm, local_bundle_width) - a
   monotonic-relationship test
   that doesn't assume linearity, using every point with a valid value
   for both variables.

2. NEAR vs FAR GROUP COMPARISON via Mann-Whitney U - the same test
   your existing fibril_spacing_combined.R analysis already uses, for
   consistency with your established methodology. Points are split into
   "near" (bottom quartile of distance-to-that-feature) and "far" (top
   quartile), and the two groups' curvature/gap distributions are
   compared.

All tomograms are pooled together for these tests (this is a per-point
relationship question, not a per-tomogram comparison), but a
per-tomogram breakdown is also printed so you can spot-check whether
any single tomogram is driving an overall result.

REQUIRES
--------
numpy, pandas, scipy
"""

import numpy as np
import pandas as pd
from scipy import stats

# ============================== CONFIG ==================================

# Reads the output of bundle_width_analysis.py, which is
# fibril_perturbation_analysis.py's output with local_bundle_width added.
# To run on the perturbation output alone, set this to
# "fibril_perturbation_analysis.csv" and remove local_bundle_width from
# METRICS below.
INPUT_CSV = "fibril_perturbation_with_bundle_width.csv"
OUTPUT_PREFIX = "perturbation_feature_correlations"

METRICS = ["curvature", "local_sustained_gap_nm", "local_bundle_width"]

QUARTILE = 0.25

# ==========================================================================


def analyze_feature_metric(df, feature_col, metric_col):
    sub = df[[feature_col, metric_col]].dropna()
    if len(sub) < 10:
        return None

    rho, p_spearman = stats.spearmanr(sub[feature_col], sub[metric_col])

    low_cut = sub[feature_col].quantile(QUARTILE)
    high_cut = sub[feature_col].quantile(1 - QUARTILE)
    near = sub[sub[feature_col] <= low_cut][metric_col]
    far = sub[sub[feature_col] >= high_cut][metric_col]

    if len(near) < 3 or len(far) < 3:
        u_stat, p_mw = np.nan, np.nan
    else:
        u_stat, p_mw = stats.mannwhitneyu(near, far, alternative="two-sided")

    return {
        "feature": feature_col.replace("dist_to_", ""),
        "metric": metric_col,
        "n_points": len(sub),
        "spearman_rho": rho,
        "spearman_p": p_spearman,
        "near_median": near.median() if len(near) > 0 else np.nan,
        "far_median": far.median() if len(far) > 0 else np.nan,
        "near_n": len(near),
        "far_n": len(far),
        "mannwhitney_p": p_mw,
    }


def main():
    df = pd.read_csv(INPUT_CSV)
    feature_cols = [c for c in df.columns if c.startswith("dist_to_")]
    metrics = [m for m in METRICS if m in df.columns]
    missing = [m for m in METRICS if m not in df.columns]
    print(f"Loaded {len(df)} points from {INPUT_CSV}")
    print(f"Feature types found: {[c.replace('dist_to_', '') for c in feature_cols]}")
    print(f"Metrics found: {metrics}")
    if missing:
        print(f"Metrics not present in this input, skipped: {missing}")
    print()

    all_results = []
    for feature_col in feature_cols:
        for metric_col in metrics:
            result = analyze_feature_metric(df, feature_col, metric_col)
            if result is not None:
                all_results.append(result)

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(f"{OUTPUT_PREFIX}_pooled.csv", index=False)

    print("=" * 100)
    print("POOLED RESULTS (all tomograms combined)")
    print("=" * 100)
    for _, row in results_df.iterrows():
        sig_spear = "*" if row["spearman_p"] < 0.05 else " "
        sig_mw = "*" if pd.notna(row["mannwhitney_p"]) and row["mannwhitney_p"] < 0.05 else " "
        print(f"\n{row['feature']} vs {row['metric']}  (n={row['n_points']})")
        print(f"  Spearman rho={row['spearman_rho']:+.3f}, p={row['spearman_p']:.4f} {sig_spear}")
        print(f"  Near (bottom {int(QUARTILE*100)}%) median={row['near_median']:.4f} "
              f"(n={row['near_n']}), Far (top {int(QUARTILE*100)}%) median={row['far_median']:.4f} "
              f"(n={row['far_n']})")
        print(f"  Mann-Whitney U p={row['mannwhitney_p']:.4f} {sig_mw}")

    print()
    print("* = p < 0.05. Spearman rho: negative means the metric DECREASES as")
    print("distance-to-feature INCREASES, i.e. the metric is elevated NEAR the")
    print("feature. Positive means the opposite (elevated far from the feature).")
    print()

    print("=" * 100)
    print("PER-TOMOGRAM BREAKDOWN (spearman rho only, for spotting tomogram-driven effects)")
    print("=" * 100)
    per_tomo_results = []
    for tomogram, tomo_df in df.groupby("tomogram"):
        for feature_col in feature_cols:
            for metric_col in metrics:
                result = analyze_feature_metric(tomo_df, feature_col, metric_col)
                if result is not None:
                    result["tomogram"] = tomogram
                    per_tomo_results.append(result)
    per_tomo_df = pd.DataFrame(per_tomo_results)
    if len(per_tomo_df) > 0:
        per_tomo_df.to_csv(f"{OUTPUT_PREFIX}_per_tomogram.csv", index=False)
        pivot = per_tomo_df.pivot_table(
            index=["feature", "metric"], columns="tomogram", values="spearman_rho"
        )
        print(pivot.to_string())

    print()
    print(f"Saved pooled results to {OUTPUT_PREFIX}_pooled.csv")
    print(f"Saved per-tomogram breakdown to {OUTPUT_PREFIX}_per_tomogram.csv")


if __name__ == "__main__":
    main()
