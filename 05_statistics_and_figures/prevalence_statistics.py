#!/usr/bin/env python3
"""
prevalence_statistics.py

Statistical testing of ultrastructural feature prevalence in annotated
tomograms, for the two prevalence tables reported in the mitochondrial
chapter.

TABLE 1 compares feature prevalence between sample types (Control,
wtTau, 2xTau). TABLE 2 compares feature prevalence between proximity
tiers within the 2xTau sample (Within an inclusion, Near an inclusion,
Distal from any inclusion).

Both are counts of tomograms in which a feature was annotated, out of
the total number of tomograms in that group, so differences are tested
by chi-squared test of independence on the underlying counts. Fisher's
exact test is used for pairwise comparisons and wherever expected cell
counts fall below five.

A separate comparison is made for granule-containing mitochondria as a
PROPORTION OF ALL MITOCHONDRIA observed, rather than as a proportion of
all tomograms. This controls for the differing overall prevalence of
mitochondria between groups, and is the comparison reported in the
text.

REQUIRES
--------
numpy, scipy
"""

import itertools
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact

# ============================== DATA ====================================
# Counts of tomograms in which each feature was annotated.

SAMPLE_TYPE_TOTALS = {"Control": 258, "wtTau": 242, "2xTau": 404}

SAMPLE_TYPE_COUNTS = {
    "Endoplasmic Reticulum":    (121, 118, 208),
    "Microtubules":             (129, 146, 171),
    "Projections":              (111, 138, 138),
    "Ribosomes":                ( 89,  89, 179),
    "Mitochondria":             ( 67,  79, 126),
    "Mitochondria w/ Granules": ( 19,  17,  68),
    "Cytoskeleton":             ( 78,  51,  91),
    "Myelin":                   ( 31,  17,  54),
    "Debris":                   (  1,  24,  39),
    "Synapse":                  ( 10,  17,  37),
    "Large Aggregate":          (  0,   0,   9),
}

PROXIMITY_TOTALS = {"Within": 31, "Near": 178, "Distal": 195}

PROXIMITY_COUNTS = {
    "Endoplasmic Reticulum":    (19, 90, 99),
    "Microtubules":             (11, 64, 96),
    "Projections":              ( 2, 78, 58),
    "Ribosomes":                (18, 65, 96),
    "Mitochondria":             ( 5, 71, 50),
    "Mitochondria w/ Granules": ( 4, 20, 44),
    "Cytoskeleton":             ( 2, 25, 64),
    "Myelin":                   ( 0, 28, 26),
    "Debris":                   ( 1, 17, 21),
    "Synapse":                  ( 1, 17, 19),
    "Large Aggregate":          ( 0,  5,  4),
}

ALPHA = 0.05

# ==========================================================================


def test_prevalence(counts, totals, table_name):
    """Chi-squared test of independence for each feature, on the counts
    of tomograms with and without that feature in each group."""
    group_names = list(totals.keys())
    group_totals = list(totals.values())

    print("=" * 78)
    print(table_name)
    print(f"Groups: " + ", ".join(f"{n} (n={t})" for n, t in zip(group_names, group_totals)))
    print("=" * 78)
    print(f"{'Feature':<28}{'chi2':>9}{'p':>14}   {'':<8}{'low expected':>14}")
    print("-" * 78)

    results = {}
    for feature, present in counts.items():
        absent = [tot - p for tot, p in zip(group_totals, present)]
        table = np.array([list(present), absent])
        chi2, p, _, expected = chi2_contingency(table)
        low = (expected < 5).any()
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < ALPHA else "ns"
        note = "yes" if low else ""
        print(f"{feature:<28}{chi2:>9.2f}{p:>14.3e}   {star:<8}{note:>14}")
        results[feature] = p

    print()
    print("Features with low expected counts should be interpreted via")
    print("Fisher's exact test on the relevant pairwise comparison rather")
    print("than the chi-squared approximation.")
    print()
    return results


def test_granule_proportion(counts, totals, group_names, table_name):
    """Granule-containing mitochondria as a proportion of ALL
    mitochondria observed, rather than of all tomograms. Controls for
    differing mitochondrial prevalence between groups."""
    with_granules = list(counts["Mitochondria w/ Granules"])
    without_granules = list(counts["Mitochondria"])

    print("=" * 78)
    print(f"{table_name}: granule-containing mitochondria as a proportion")
    print("of all mitochondria observed")
    print("=" * 78)

    for name, wg, wo in zip(group_names, with_granules, without_granules):
        total = wg + wo
        pct = 100 * wg / total if total else float("nan")
        lo, hi = wilson_interval(wg, total)
        print(f"  {name:<10} {wg:>3} / {total:<4} = {pct:5.1f}%   "
              f"95% CI {lo:5.1f} to {hi:5.1f}%")

    table = np.array([with_granules, without_granules])
    chi2, p, _, expected = chi2_contingency(table)
    print()
    print(f"  Overall chi-squared: chi2 = {chi2:.2f}, p = {p:.4f}")
    if (expected < 5).any():
        print("  (low expected counts present, interpret with caution)")

    print()
    print("  Pairwise, Fisher's exact test with Bonferroni correction:")
    n_comparisons = len(list(itertools.combinations(range(len(group_names)), 2)))
    for i, j in itertools.combinations(range(len(group_names)), 2):
        _, p_pair = fisher_exact([[with_granules[i], without_granules[i]],
                                   [with_granules[j], without_granules[j]]])
        p_corrected = min(1.0, p_pair * n_comparisons)
        star = "*" if p_corrected < ALPHA else "ns"
        print(f"    {group_names[i]:<10} vs {group_names[j]:<10} "
              f"p = {p_pair:.4f}, corrected p = {p_corrected:.4f}  {star}")
    print()


def wilson_interval(successes, total, z=1.96):
    """Wilson score interval for a binomial proportion, expressed as a
    percentage. More reliable than the normal approximation at small n."""
    if total == 0:
        return float("nan"), float("nan")
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    spread = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return 100 * (centre - spread), 100 * (centre + spread)


def main():
    test_prevalence(SAMPLE_TYPE_COUNTS, SAMPLE_TYPE_TOTALS,
                    "TABLE 1: feature prevalence by sample type")

    test_granule_proportion(SAMPLE_TYPE_COUNTS, SAMPLE_TYPE_TOTALS,
                            list(SAMPLE_TYPE_TOTALS.keys()),
                            "Sample type")

    test_prevalence(PROXIMITY_COUNTS, PROXIMITY_TOTALS,
                    "TABLE 2: feature prevalence by proximity to an inclusion (2xTau only)")

    test_granule_proportion(PROXIMITY_COUNTS, PROXIMITY_TOTALS,
                            list(PROXIMITY_TOTALS.keys()),
                            "Proximity tier")


if __name__ == "__main__":
    main()
