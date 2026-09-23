#!/usr/bin/env python3
"""
bundle_width_analysis.py

Extends fibril_perturbation_analysis.py with a local bundle-width metric.

Where fibril_perturbation_analysis.py measures the gap to the single
nearest qualifying neighbour, this asks a different question: how many
parallel fibrils are travelling alongside this point at all? A fibril
running alone through the volume and a fibril embedded in the middle of
a tightly packed bundle can have identical nearest-neighbour gaps, but
very different local environments. Bundle width distinguishes them.

DEFINITION
----------
For each point on each traced fibril, a neighbouring fibril is counted
towards that point's bundle width if all three of the following hold:

  1. Its local axis lies within ANGLE_THRESH degrees of this point's own
     local axis (same parallelism criterion as the sustained-gap
     measure).
  2. Its nearest point lies within MAX_SEARCH voxels of this point.
  3. That nearest point is offset along the neighbour's own axis by no
     more than WINDOW_VOXELS, so a fibril that merely crosses this one
     at a shallow angle is not counted as running alongside it.

Local axes are estimated at each point by principal component analysis
over an arc-length window of WINDOW_VOXELS on either side, the same
approach used throughout the fibril analyses.

INPUT
-----
fibril_perturbation_analysis.csv, as produced by
fibril_perturbation_analysis.py (requires columns: tomogram, fibril_id,
point_index, x, y, z).

OUTPUT
------
fibril_perturbation_with_bundle_width.csv - the input with an added
local_bundle_width column, ready to correlate against distance to
ultrastructural features in correlate_perturbation_with_features.py.

REQUIRES
--------
numpy, pandas, scipy
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import stats

# ============================== CONFIG ==================================

INPUT_CSV = "fibril_perturbation_analysis.csv"
OUTPUT_CSV = "fibril_perturbation_with_bundle_width.csv"

PIXEL_SIZE_ANGSTROM = 9.92      # 2.48 A/pixel at bin 4

# Arc-length window for local axis estimation, and the maximum
# longitudinal offset permitted for a neighbour to count as running
# alongside rather than crossing. 20 nm.
WINDOW_NM = 20.0
WINDOW_VOXELS = (WINDOW_NM * 10) / PIXEL_SIZE_ANGSTROM

# Maximum centreline separation, in voxels, at which a neighbour is
# considered. 100 voxels at bin 4 is approximately 99 nm, which spans
# roughly three shells of neighbours at the observed 31 nm sustained
# spacing.
MAX_SEARCH = 100.0

# Parallelism criterion, matching the sustained-gap measure.
ANGLE_THRESH = 8.0

# ==========================================================================


def local_axis(points, cumulative, idx, window_voxels):
    """Estimates the local trace direction at points[idx] by PCA over a
    window of points within window_voxels arc-length either side.
    Returns (unit_axis, window_points), or None if too few points fall
    within the window."""
    center = cumulative[idx]
    mask = np.abs(cumulative - center) <= window_voxels
    window_points = points[mask]
    if len(window_points) < 2:
        return None
    centered = window_points - window_points.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axis = vt[0]
    norm = np.linalg.norm(axis)
    if norm == 0:
        return None
    return axis / norm, window_points


def build_fibril_data(tomogram_df):
    """Groups a tomogram's points by fibril, computes cumulative
    arc-length along each, and builds a KD-tree per fibril."""
    fibril_data = {}
    for fid, group in tomogram_df.groupby("fibril_id"):
        points = group.sort_values("point_index")[["x", "y", "z"]].to_numpy()
        if len(points) < 2:
            continue
        seg_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
        cumulative = np.concatenate([[0], np.cumsum(seg_lengths)])
        fibril_data[fid] = {
            "points": points,
            "cum": cumulative,
            "tree": cKDTree(points),
        }
    return fibril_data


def compute_bundle_widths(fibril_data):
    """Returns {fibril_id: list of bundle widths, one per point}."""
    widths_by_fibril = {}

    for fid, data in fibril_data.items():
        points, cumulative = data["points"], data["cum"]
        widths = []

        for i in range(len(points)):
            own = local_axis(points, cumulative, i, WINDOW_VOXELS)
            if own is None:
                widths.append(np.nan)
                continue
            own_axis, _ = own

            count = 0
            for gid, gdata in fibril_data.items():
                if gid == fid:
                    continue

                dist, nearest_idx = gdata["tree"].query(points[i])
                if dist > MAX_SEARCH:
                    continue

                neighbour = local_axis(gdata["points"], gdata["cum"],
                                       nearest_idx, WINDOW_VOXELS)
                if neighbour is None:
                    continue
                neighbour_axis, neighbour_window = neighbour
                if len(neighbour_window) < 3:
                    continue

                cos_angle = abs(float(np.dot(own_axis, neighbour_axis)))
                angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                if angle > ANGLE_THRESH:
                    continue

                offset = points[i] - gdata["points"][nearest_idx]
                along = float(np.dot(offset, neighbour_axis))
                if abs(along) > WINDOW_VOXELS:
                    continue

                count += 1

            widths.append(count)

        widths_by_fibril[fid] = widths

    return widths_by_fibril


def main():
    df = pd.read_csv(INPUT_CSV)
    print(f"Loaded {len(df)} points from {INPUT_CSV}")
    print(f"Parallelism threshold: {ANGLE_THRESH} degrees")
    print(f"Maximum centreline separation: {MAX_SEARCH} voxels "
          f"({MAX_SEARCH * PIXEL_SIZE_ANGSTROM / 10:.0f} nm)")
    print(f"Axis window and longitudinal offset limit: {WINDOW_NM} nm "
          f"({WINDOW_VOXELS:.1f} voxels)")
    print()

    all_rows = []
    for tomogram, tomogram_df in df.groupby("tomogram"):
        fibril_data = build_fibril_data(tomogram_df)
        widths_by_fibril = compute_bundle_widths(fibril_data)

        for fid, widths in widths_by_fibril.items():
            for i, w in enumerate(widths):
                all_rows.append({
                    "tomogram": tomogram,
                    "fibril_id": fid,
                    "point_index": i,
                    "local_bundle_width": w,
                })
        print(f"  {tomogram}: {len(fibril_data)} fibrils processed")

    width_df = pd.DataFrame(all_rows)
    merged = df.merge(width_df, on=["tomogram", "fibril_id", "point_index"],
                      how="left")
    merged.to_csv(OUTPUT_CSV, index=False)

    print()
    print("Local bundle width distribution:")
    print(merged["local_bundle_width"].describe().to_string())
    print()
    print("Counts by value:")
    print(merged["local_bundle_width"].value_counts().sort_index().to_string())

    # relationship between curvature and local packing, independent of
    # any ultrastructural feature
    sub = merged[["curvature", "local_sustained_gap_nm", "tomogram"]].dropna()
    if len(sub) > 10:
        rho, p = stats.spearmanr(sub["curvature"], sub["local_sustained_gap_nm"])
        print()
        print(f"Curvature vs local sustained gap, pooled (n={len(sub)}): "
              f"rho={rho:+.3f}, p={p:.6f}")
        for tomo, tdf in sub.groupby("tomogram"):
            r, pv = stats.spearmanr(tdf["curvature"], tdf["local_sustained_gap_nm"])
            print(f"  {tomo}: rho={r:+.3f}, p={pv:.4f}, n={len(tdf)}")

    print()
    print(f"Wrote {OUTPUT_CSV}")
    print("Add 'local_bundle_width' to METRICS in "
          "correlate_perturbation_with_features.py to test it against "
          "distance to each ultrastructural feature.")


if __name__ == "__main__":
    main()
