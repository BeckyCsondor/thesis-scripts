#!/usr/bin/env python3
"""
fibril_perturbation_analysis.py

For every point along every traced fibril, computes:
  - distance to the nearest ultrastructural feature of each type
    (ER, mito: nearest boundary point. ribos: distance to that
    ribosome's center of mass. ribocluster: nearest boundary point.
    microtubule: nearest confident microtubule pick, treated the same
    way as a boundary-type feature)
  - local curvature at that point (from the trace's own geometry)
  - LOCAL SUSTAINED GAP to the nearest genuinely-parallel neighboring
    fibril, using your own validated definition (matching your existing
    fibril_spacing_combined.R methodology): pairs within
    ANGLE_THRESH_DEG of parallel, gap measured as perpendicular
    distance minus FIBRIL_WIDTH_NM (center-to-center -> surface-to-
    surface). The difference here is that this is made LOCAL - instead
    of one rigid whole-fibril axis and a single global overlap check,
    each point's own local axis is estimated from a small WINDOW_NM
    window of nearby trace points, and compared against the nearest
    OTHER fibril's own local axis in the same way - so the gap value
    can change along a fibril's length, and can be correlated against
    curvature and feature proximity at each point.

Output is one row per fibril point: curvature, local_sustained_gap_nm,
and a distance-to-nearest column for each feature type present.

FEATURE FILE NAMING CONVENTION - ULTRASTRUCTURAL FEATURES
------------------------------------------------------------
Flat files directly in ULTRASTRUCTURE_DIR, named:
    {TomogramName}_{FeatureName}_points.txt
Feature TYPE is inferred from the feature name: anything starting with
"ribos" (but not "ribocluster") is center-of-mass distance; everything
else is nearest-boundary-point distance.

MICROTUBULE FEATURE - separate config, since these come from your
pytme confident-microtubule .tsv files, not the ultrastructure point
files
------------------------------------------------------------------
MICROTUBULE_TSV_FILES: one entry per tomogram, pointing at that
tomogram's confident microtubule .tsv (from remove_microtubule_overlap.py's
output). Treated as a single nearest-boundary-type feature ("microtubule").

FIBRIL TRACE INPUT
--------------------
FIBRIL_TRACE_FILES: one entry per tomogram, "object contour x y z"
format - use your real guided traces (assign_picks_to_manual_guides.py
output), or a spline_trace_and_pick.py dense output for smoother
curvature/gap values.

REQUIRES
--------
numpy, pandas, scipy
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import glob
import os
import re

# ============================== CONFIG ==================================

ULTRASTRUCTURE_DIR = "tomograms/ultrastructural_features"

FIBRIL_TRACE_FILES = {
    # "Position018": "polarity_withguides/Position018_guided_trace_points.txt",
}

MICROTUBULE_TSV_FILES = {
    # "Position018": "coordinates2/Position018_microtubule_confident.tsv",
}

PIXEL_SIZE_ANGSTROM = 9.92

# Local sustained-gap parameters - matching your validated
# fibril_spacing_combined.R definition, made local via WINDOW_NM.
WINDOW_NM = 20.0
ANGLE_THRESH_DEG = 8.0
FIBRIL_WIDTH_NM = 6.0
MAX_NEIGHBOR_SEARCH_VOXELS = 100.0

OUTPUT_PREFIX = "fibril_perturbation_analysis"

# ==========================================================================


def discover_feature_files(ultrastructure_dir, tomogram_name):
    """Finds every {tomogram}_{feature}_points.txt file for one
    tomogram, and classifies each by measurement type."""
    pattern = os.path.join(ultrastructure_dir, f"{tomogram_name}_*_points.txt")
    files = glob.glob(pattern)

    features = []
    for f in files:
        base = os.path.basename(f)
        m = re.match(rf"^{re.escape(tomogram_name)}_(.+)_points\.txt$", base)
        if not m:
            continue
        feature_instance_name = m.group(1)  # e.g. "ER1", "mito2", "ribos", "ribocluster"

        if feature_instance_name.lower().startswith("ribos") and \
           "cluster" not in feature_instance_name.lower():
            measurement = "center_of_mass"
        else:
            measurement = "nearest_boundary"

        # strip trailing digits to get the general feature TYPE
        # (ER1, ER2, ER3 -> all type "ER")
        feature_type = re.sub(r"\d+$", "", feature_instance_name)

        features.append({
            "instance_name": feature_instance_name,
            "type": feature_type,
            "measurement": measurement,
            "path": f,
        })
    return features


def load_feature_points(path):
    cols = ["object", "contour", "x", "y", "z"]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)
    return df[["x", "y", "z"]].to_numpy()


def load_fibril_traces(path):
    cols = ["object", "contour", "x", "y", "z"]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)
    df["fibril_id"] = df["object"].astype(str) + "_" + df["contour"].astype(str)
    return df


def local_axis_window(points, cumulative, idx, window_voxels):
    """Estimates the local trace direction at points[idx], using a
    window of nearby points (within window_voxels arc-length distance
    on either side), via PCA/SVD - same principle as every other local
    direction estimate tonight, just parameterized in real physical
    units (window_nm) here. Returns (unit_axis, window_points) or None
    if too few points fall in the window."""
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


def compute_local_sustained_gaps(fibril_points, pixel_size_angstrom, window_nm,
                                  angle_thresh_deg, fibril_width_nm,
                                  max_search_voxels):
    """For every point on every fibril, finds the nearest OTHER fibril
    that is locally parallel (within angle_thresh_deg) with enough
    local extent to count as genuinely running alongside (not just
    touching at one point), and returns the perpendicular gap to it
    (minus fibril_width_nm) in nm - your validated sustained-gap
    definition, made local. NaN where no qualifying neighbor exists at
    that point. Returns {fibril_id: gap_array_nm}."""
    window_voxels = (window_nm * 10) / pixel_size_angstrom
    width_voxels = (fibril_width_nm * 10) / pixel_size_angstrom

    fibril_data = {}
    for fid, points in fibril_points.items():
        if len(points) < 2:
            fibril_data[fid] = {"points": points, "cumulative": np.zeros(len(points)), "tree": None}
            continue
        seg_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
        cumulative = np.concatenate([[0], np.cumsum(seg_lengths)])
        fibril_data[fid] = {"points": points, "cumulative": cumulative, "tree": cKDTree(points)}

    fids = list(fibril_points.keys())
    results = {}

    for fid in fids:
        data_f = fibril_data[fid]
        points_f, cumulative_f = data_f["points"], data_f["cumulative"]
        n = len(points_f)
        gaps = np.full(n, np.nan)

        for i in range(n):
            axis_result_f = local_axis_window(points_f, cumulative_f, i, window_voxels)
            if axis_result_f is None:
                continue
            axis_f, _ = axis_result_f

            best_gap_nm = None
            for gid in fids:
                if gid == fid:
                    continue
                data_g = fibril_data[gid]
                if data_g["tree"] is None:
                    continue

                dist, nearest_idx = data_g["tree"].query(points_f[i])
                if dist > max_search_voxels:
                    continue

                axis_result_g = local_axis_window(data_g["points"], data_g["cumulative"],
                                                    nearest_idx, window_voxels)
                if axis_result_g is None:
                    continue
                axis_g, window_points_g = axis_result_g

                if len(window_points_g) < 3:
                    continue  # not enough local extent to count as "running alongside"

                cos_angle = abs(float(np.dot(axis_f, axis_g)))
                angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
                if angle_deg > angle_thresh_deg:
                    continue

                nearest_g_point = data_g["points"][nearest_idx]
                v = points_f[i] - nearest_g_point
                along_component = float(np.dot(v, axis_g))
                # reject if point i sits genuinely PAST the end of G's local
                # window, not beside it - this is what "overlap" meant in the
                # original whole-fibril definition, applied locally: a large
                # along-axis offset from G's nearest point means G doesn't
                # actually extend alongside point i, even if the raw nearest-
                # point distance looked small
                if abs(along_component) > window_voxels:
                    continue
                perp_dist_voxels = np.linalg.norm(v - along_component * axis_g)

                gap_voxels = perp_dist_voxels - width_voxels
                gap_nm = (gap_voxels * pixel_size_angstrom) / 10.0

                if best_gap_nm is None or gap_nm < best_gap_nm:
                    best_gap_nm = gap_nm

            gaps[i] = best_gap_nm if best_gap_nm is not None else np.nan

        results[fid] = gaps

    return results



    """Three-point discrete curvature estimate at each interior point.
    Returns an array the same length as points, NaN at both endpoints."""
def compute_curvature(points):
    """Three-point discrete curvature estimate at each interior point.
    Returns an array the same length as points, NaN at both endpoints."""
    n = len(points)
    curvature = np.full(n, np.nan)
    for i in range(1, n - 1):
        p0, p1, p2 = points[i - 1], points[i], points[i + 1]
        a = np.linalg.norm(p1 - p0)
        b = np.linalg.norm(p2 - p1)
        c = np.linalg.norm(p2 - p0)
        if a == 0 or b == 0 or c == 0:
            continue
        # area of the triangle via cross product, then standard
        # curvature-from-circumradius formula: curvature = 4*Area / (a*b*c)
        cross = np.cross(p1 - p0, p2 - p0)
        area = 0.5 * np.linalg.norm(cross)
        if area == 0:
            curvature[i] = 0.0  # perfectly straight, three collinear points
        else:
            curvature[i] = (4 * area) / (a * b * c)
    return curvature


def load_microtubule_points(path):
    df = pd.read_csv(path, sep="\t")
    return df[["x", "y", "z"]].to_numpy()


def main():
    if not FIBRIL_TRACE_FILES:
        print("FIBRIL_TRACE_FILES is empty - add an entry per tomogram, then rerun.")
        print()
        print("Showing what feature files WOULD be picked up per tomogram,")
        print("so you can sanity-check the naming/classification now:")
        print()
        tomogram_names = set()
        for f in glob.glob(os.path.join(ULTRASTRUCTURE_DIR, "*_points.txt")):
            base = os.path.basename(f)
            m = re.match(r"^(Position\d+)_.+_points\.txt$", base)
            if m:
                tomogram_names.add(m.group(1))
        for name in sorted(tomogram_names):
            feats = discover_feature_files(ULTRASTRUCTURE_DIR, name)
            print(f"{name}:")
            for feat in feats:
                print(f"  {feat['instance_name']} (type={feat['type']}, "
                      f"measurement={feat['measurement']}) <- {feat['path']}")
            if name in MICROTUBULE_TSV_FILES:
                print(f"  microtubule (type=microtubule, measurement=nearest_boundary) "
                      f"<- {MICROTUBULE_TSV_FILES[name]}")
        return

    all_results = []

    for tomogram_name, trace_path in FIBRIL_TRACE_FILES.items():
        print(f"Processing {tomogram_name}...")
        fibrils = load_fibril_traces(trace_path)
        feature_specs = discover_feature_files(ULTRASTRUCTURE_DIR, tomogram_name)

        feature_lookup = {}
        for feat in feature_specs:
            points = load_feature_points(feat["path"])
            if feat["measurement"] == "center_of_mass":
                reference = points.mean(axis=0)
            else:
                reference = cKDTree(points)
            feature_lookup.setdefault(feat["type"], []).append({
                "instance_name": feat["instance_name"],
                "measurement": feat["measurement"],
                "reference": reference,
            })

        if tomogram_name in MICROTUBULE_TSV_FILES:
            mt_points = load_microtubule_points(MICROTUBULE_TSV_FILES[tomogram_name])
            if len(mt_points) > 0:
                feature_lookup.setdefault("microtubule", []).append({
                    "instance_name": "microtubule",
                    "measurement": "nearest_boundary",
                    "reference": cKDTree(mt_points),
                })

        feature_types = sorted(feature_lookup.keys())
        print(f"  Feature types found: {feature_types}")

        # ---- local sustained-gap between parallel fibrils ----
        fibril_points_dict = {
            fid: group[["x", "y", "z"]].to_numpy()
            for fid, group in fibrils.groupby("fibril_id", sort=False)
        }
        print(f"  Computing local sustained gaps across "
              f"{len(fibril_points_dict)} fibrils (this can take a while)...")
        gap_results = compute_local_sustained_gaps(
            fibril_points_dict, PIXEL_SIZE_ANGSTROM, WINDOW_NM,
            ANGLE_THRESH_DEG, FIBRIL_WIDTH_NM, MAX_NEIGHBOR_SEARCH_VOXELS
        )

        for fibril_id, points in fibril_points_dict.items():
            curvature = compute_curvature(points)
            local_gap_nm = gap_results[fibril_id]

            row_data = {
                "tomogram": tomogram_name,
                "fibril_id": fibril_id,
                "point_index": np.arange(len(points)),
                "x": points[:, 0], "y": points[:, 1], "z": points[:, 2],
                "curvature": curvature,
                "local_sustained_gap_nm": local_gap_nm,
            }

            for ftype in feature_types:
                best_dist = np.full(len(points), np.inf)
                for instance in feature_lookup[ftype]:
                    if instance["measurement"] == "center_of_mass":
                        d = np.linalg.norm(points - instance["reference"], axis=1)
                    else:
                        d, _ = instance["reference"].query(points)
                    best_dist = np.minimum(best_dist, d)
                row_data[f"dist_to_{ftype}"] = best_dist

            all_results.append(pd.DataFrame(row_data))

        print(f"  {len(fibril_points_dict)} fibrils processed")

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(f"{OUTPUT_PREFIX}.csv", index=False)

    print()
    print(f"Wrote {len(combined)} fibril points across "
          f"{len(FIBRIL_TRACE_FILES)} tomogram(s) to {OUTPUT_PREFIX}.csv")
    print()
    n_with_gap = combined["local_sustained_gap_nm"].notna().sum()
    print(f"{n_with_gap} / {len(combined)} points had a qualifying sustained-parallel "
          f"neighbor within {MAX_NEIGHBOR_SEARCH_VOXELS} voxels "
          f"(angle <= {ANGLE_THRESH_DEG} deg)")
    print()
    print("Each row: one fibril point, curvature, local_sustained_gap_nm (to the")
    print("nearest genuinely-parallel neighboring fibril, or NaN if none nearby),")
    print("and distance to the nearest instance of each feature type.")
    print("Bin/plot curvature and local_sustained_gap_nm against dist_to_<feature>")
    print("columns to look for association with proximity to ultrastructure.")


if __name__ == "__main__":
    main()
