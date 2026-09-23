#!/usr/bin/env python3
"""
polarity_and_bundle_analysis.py

Combines polarity_analysis.py and bundle_clustering.py into a single
pipeline. Two real benefits beyond convenience:

1. bundle_clustering.py never had a "don't bundle across different
   tomograms" restriction (unlike polarity_analysis.py's pairwise
   comparison, which was specifically fixed for this earlier) -
   comparing spatial bundles across unrelated tomograms would be just
   as meaningless as comparing polarity across them. Since bundling
   now happens INSIDE the same per-tomogram loop as everything else,
   this is naturally enforced, not something that could be forgotten.

2. Filament centroids for bundling are computed directly from each
   filament's ASSIGNED PICKS in memory, rather than by re-loading and
   re-matching a saved CSV's filament_id strings against a contour
   file a second time - which is exactly the "TomogramName::" prefix
   mismatch bug that bit both make_polarity_model.py and
   bundle_clustering.py tonight. Keeping everything in one continuous
   in-memory pipeline avoids that whole bug category going forward.

WHAT IT DOES, PER TOMOGRAM
-----------------------------
1. Loads picks + traced contours, assigns picks to filaments, computes
   each filament's direction/internal_agreement (same as
   polarity_analysis.py).
2. Computes pairwise same/opposite/ambiguous polarity for that
   tomogram's filaments only.
3. Computes spatial/parallelism bundles for that tomogram's filaments
   only, and reports same/opposite/mixed polarity per bundle.
Then combines all tomograms' results for overall summary tallies.

REQUIRES
--------
numpy, pandas, scipy
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

# ============================== CONFIG ==================================

TOMOGRAMS = [
    {
        "name": "Position018",
        "tsv": "coordinates2/Position018_bin4_2_no_mt_overlap.tsv",
        "contours": "polarity_withguides/Position018_guided_trace_points.txt",
    },
    # Add Position019, Position022, Position023 here in the same format
]

OUTPUT_PREFIX = "polarity_withguides/combined_analysis"

EULER_CONVENTION = "zyx"
REFERENCE_AXIS = np.array([0.0, 1.0, 0.0])

MAX_ASSIGNMENT_DISTANCE = 15.0

MIN_PICKS_PER_FILAMENT = 5
MIN_AGREEMENT = 0.6

MAX_BUNDLE_DISTANCE = 100.0
MIN_PARALLELISM_INDEX = 0.95

# ==========================================================================


def load_picks(tsv_path):
    return pd.read_csv(tsv_path, sep="\t")


def load_contours(path):
    cols = ["object", "contour", "x", "y", "z"]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)
    df["filament_id"] = df["object"].astype(str) + "_" + df["contour"].astype(str)
    return df


def sanity_check(picks, convention, reference_axis, tomogram_name):
    print(f"[{tomogram_name}] Sanity check - simple single-axis rotations:")
    subset = picks[
        (np.isclose(picks["euler_z"], 0, atol=0.5) & np.isclose(picks["euler_y"], 0, atol=0.5))
        | (np.isclose(picks["euler_z"], 0, atol=0.5) & np.isclose(picks["euler_x"], 0, atol=0.5))
        | (np.isclose(picks["euler_y"], 0, atol=0.5) & np.isclose(picks["euler_x"], 0, atol=0.5))
    ].head(4)
    if subset.empty:
        print("  (no simple single-axis rows found in top of file)")
        return
    angles = subset[["euler_z", "euler_y", "euler_x"]].to_numpy()
    rot = Rotation.from_euler(convention.upper(), angles, degrees=True)
    vecs = rot.apply(reference_axis)
    for (_, row), vec in zip(subset.iterrows(), vecs):
        print(f"  ({row['x']:.0f},{row['y']:.0f},{row['z']:.0f}) "
              f"euler=({row['euler_z']:.1f},{row['euler_y']:.1f},{row['euler_x']:.1f}) "
              f"-> dir=({vec[0]:+.2f},{vec[1]:+.2f},{vec[2]:+.2f})")
    print()


def euler_angles_to_vectors(picks, convention, reference_axis):
    angles = picks[["euler_z", "euler_y", "euler_x"]].to_numpy()
    rot = Rotation.from_euler(convention.upper(), angles, degrees=True)
    vectors = rot.apply(reference_axis)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return vectors / norms


def assign_picks_to_filaments(picks, contours, max_dist):
    tree = cKDTree(contours[["x", "y", "z"]].to_numpy())
    pick_xyz = picks[["x", "y", "z"]].to_numpy()
    dist, idx = tree.query(pick_xyz)
    filament_id = contours["filament_id"].to_numpy()[idx]
    picks = picks.copy()
    picks["nearest_filament"] = filament_id
    picks["assignment_distance"] = dist
    picks.loc[picks["assignment_distance"] > max_dist, "nearest_filament"] = None
    return picks


def cluster_filament_direction(vectors, n_iter=5):
    if len(vectors) == 0:
        return None, 0.0
    mean_vec = vectors[0].copy()
    for _ in range(n_iter):
        signs = np.sign(vectors @ mean_vec)
        signs[signs == 0] = 1.0
        aligned = vectors * signs[:, None]
        mean_vec = aligned.mean(axis=0)
        norm = np.linalg.norm(mean_vec)
        if norm > 0:
            mean_vec = mean_vec / norm
    agreement = np.mean((vectors @ mean_vec) > 0)
    agreement = max(agreement, 1 - agreement)
    return mean_vec, agreement


def compute_pairwise(filament_df):
    ids = filament_df["filament_id"].tolist()
    vecs = filament_df[["mean_dir_x", "mean_dir_y", "mean_dir_z"]].to_numpy()
    rows = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            dot = float(np.clip(vecs[i] @ vecs[j], -1.0, 1.0))
            angle = np.degrees(np.arccos(dot))
            if angle < 45:
                label = "same polarity"
            elif angle > 135:
                label = "opposite polarity"
            else:
                label = "ambiguous (not clearly parallel or antiparallel)"
            rows.append({"filament_1": ids[i], "filament_2": ids[j],
                          "angle_deg": angle, "label": label})
    return pd.DataFrame(rows)


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, i):
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, i, j):
        ri, rj = self.find(i), self.find(j)
        if ri != rj:
            self.parent[ri] = rj


def cluster_bundles(filament_df, max_dist, min_parallelism):
    """Bundles filaments WITHIN this single tomogram's filament_df only -
    caller is responsible for only passing one tomogram's filaments at
    a time, which main() does naturally by looping per tomogram."""
    n = len(filament_df)
    centroids = filament_df[["centroid_x", "centroid_y", "centroid_z"]].to_numpy()
    directions = filament_df[["mean_dir_x", "mean_dir_y", "mean_dir_z"]].to_numpy()

    uf = UnionFind(n)
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.linalg.norm(centroids[i] - centroids[j])
            if dist > max_dist:
                continue
            dot = float(np.clip(directions[i] @ directions[j], -1.0, 1.0))
            if abs(dot) >= min_parallelism:
                uf.union(i, j)

    bundle_ids = [uf.find(i) for i in range(n)]
    bundle_sizes = pd.Series(bundle_ids).value_counts()
    remap = {old: new for new, old in enumerate(bundle_sizes.index, start=1)}
    return [remap[b] for b in bundle_ids]


def main():
    all_filament_rows = []
    all_pair_dfs = []
    all_conf_pair_dfs = []
    all_bundle_summaries = []

    for tomo in TOMOGRAMS:
        name = tomo["name"]
        print("=" * 70)
        print(f"TOMOGRAM: {name}")
        print("=" * 70)

        picks = load_picks(tomo["tsv"])
        contours = load_contours(tomo["contours"])
        print(f"Loaded {len(picks)} picks, {contours['filament_id'].nunique()} traced filaments")

        sanity_check(picks, EULER_CONVENTION, REFERENCE_AXIS, name)

        picks = assign_picks_to_filaments(picks, contours, MAX_ASSIGNMENT_DISTANCE)
        assigned = picks.dropna(subset=["nearest_filament"]).copy()
        n_unassigned = len(picks) - len(assigned)
        print(f"{len(assigned)} / {len(picks)} picks assigned to a filament "
              f"({n_unassigned} unassigned)")

        vectors = euler_angles_to_vectors(assigned, EULER_CONVENTION, REFERENCE_AXIS)
        assigned[["dir_x", "dir_y", "dir_z"]] = vectors

        results = []
        for fid, group in assigned.groupby("nearest_filament"):
            vecs = group[["dir_x", "dir_y", "dir_z"]].to_numpy()
            mean_vec, agreement = cluster_filament_direction(vecs)
            centroid = group[["x", "y", "z"]].mean().to_numpy()
            results.append({
                "tomogram": name,
                "filament_id": f"{name}::{fid}",
                "n_picks": len(group),
                "mean_dir_x": mean_vec[0], "mean_dir_y": mean_vec[1], "mean_dir_z": mean_vec[2],
                "centroid_x": centroid[0], "centroid_y": centroid[1], "centroid_z": centroid[2],
                "internal_agreement": agreement,
                "mean_score": group["score"].mean(),
            })
        filament_df = pd.DataFrame(results)
        filament_df["confident"] = (
            (filament_df["n_picks"] >= MIN_PICKS_PER_FILAMENT)
            & (filament_df["internal_agreement"] >= MIN_AGREEMENT)
        )
        all_filament_rows.append(filament_df)
        print(f"Computed direction for {len(filament_df)} filaments, "
              f"{filament_df['confident'].sum()} meet the confidence bar "
              f"(n_picks >= {MIN_PICKS_PER_FILAMENT}, agreement >= {MIN_AGREEMENT})")

        pair_df = compute_pairwise(filament_df)
        pair_df["tomogram"] = name
        all_pair_dfs.append(pair_df)

        confident_only = filament_df[filament_df["confident"]].reset_index(drop=True)
        conf_pair_df = compute_pairwise(confident_only)
        conf_pair_df["tomogram"] = name
        all_conf_pair_dfs.append(conf_pair_df)
        if len(conf_pair_df) > 0:
            print(f"High-confidence pairs: {conf_pair_df['label'].value_counts().to_dict()}")

        bundle_ids = cluster_bundles(filament_df, MAX_BUNDLE_DISTANCE, MIN_PARALLELISM_INDEX)
        filament_df["bundle_id"] = bundle_ids

        bundle_rows = []
        for bid, group in filament_df.groupby("bundle_id"):
            conf_group = group[group["confident"]]
            if len(conf_group) < 2:
                bundle_rows.append({"tomogram": name, "bundle_id": bid,
                                     "n_filaments_total": len(group),
                                     "n_confident": len(conf_group),
                                     "n_same_way": np.nan, "n_opposite_way": np.nan,
                                     "status": "not enough confident members"})
                continue
            vecs = conf_group[["mean_dir_x", "mean_dir_y", "mean_dir_z"]].to_numpy()
            ref = vecs[0]
            signs = np.sign(vecs @ ref)
            signs[signs == 0] = 1.0
            n_same, n_opp = int(np.sum(signs > 0)), int(np.sum(signs < 0))
            status = "uniform" if (n_opp == 0 or n_same == 0) else "MIXED"
            bundle_rows.append({"tomogram": name, "bundle_id": bid,
                                 "n_filaments_total": len(group), "n_confident": len(conf_group),
                                 "n_same_way": n_same, "n_opposite_way": n_opp, "status": status})
        bundle_summary = pd.DataFrame(bundle_rows)
        all_bundle_summaries.append(bundle_summary)

        n_multi = (bundle_summary["n_filaments_total"] > 1).sum()
        n_verdict = bundle_summary["status"].isin(["uniform", "MIXED"]).sum()
        print(f"Bundles: {n_multi} contain >1 filament, {n_verdict} have a polarity verdict")
        print()

    combined_filaments = pd.concat(all_filament_rows, ignore_index=True)
    combined_pairs = pd.concat(all_pair_dfs, ignore_index=True)
    combined_conf_pairs = pd.concat(all_conf_pair_dfs, ignore_index=True)
    combined_bundles = pd.concat(all_bundle_summaries, ignore_index=True)

    combined_filaments.to_csv(f"{OUTPUT_PREFIX}_filaments.csv", index=False)
    combined_pairs.to_csv(f"{OUTPUT_PREFIX}_pairwise.csv", index=False)
    combined_conf_pairs.to_csv(f"{OUTPUT_PREFIX}_pairwise_highconfidence.csv", index=False)
    combined_bundles.to_csv(f"{OUTPUT_PREFIX}_bundles.csv", index=False)

    print("=" * 70)
    print(f"COMBINED RESULTS ACROSS {len(TOMOGRAMS)} TOMOGRAM(S)")
    print("(all comparisons/bundling stayed within-tomogram individually -")
    print(" these are just the combined tallies)")
    print("=" * 70)
    print()
    print("All filament pairs, combined:")
    print(combined_pairs["label"].value_counts().to_string())
    print()
    print("High-confidence filament pairs, combined:")
    if len(combined_conf_pairs) > 0:
        print(combined_conf_pairs["label"].value_counts().to_string())
    print()
    print("Bundle polarity verdicts, combined:")
    verdict_bundles = combined_bundles[combined_bundles["status"].isin(["uniform", "MIXED"])]
    if len(verdict_bundles) > 0:
        print(verdict_bundles["status"].value_counts().to_string())
    else:
        print("No bundles had enough confident members for a verdict.")
    print()
    print(f"Saved combined outputs with prefix {OUTPUT_PREFIX}")


if __name__ == "__main__":
    main()
