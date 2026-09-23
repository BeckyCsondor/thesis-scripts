#!/usr/bin/env python3
"""
assign_picks_to_manual_guides.py

Combines the best of both worlds: you manually trace a SPARSE guide path
along each fibril in 3dmod (just enough clicks to mark which visual
streak is which fibril - precision doesn't matter much here, since a
human eye trivially tells bundled fibrils apart where the automated
tracers tonight kept failing). This script then does the precise work:
it snaps your real, high-quality pytme picks onto whichever guide
contour they're nearest to, and orders them correctly along that
contour's length - producing a clean, properly-ordered trace built from
your actual pick positions (not the rough manual clicks themselves).

WORKFLOW
--------
1. In 3dmod, trace ONE open contour per fibril, sparsely (e.g. a click
   every 20-40 voxels is plenty) - just enough to unambiguously mark
   the path, especially through any tight bundles. Save the model.
2. Export it:
       module load imod
       model2point -object -float your_guide_model.mod guide_points.txt
3. Set GUIDE_POINTS_FILE and TSV_PATH below, run this script.
4. Feed the output into spline_trace_and_pick.py, same as any other
   trace output tonight.

HOW ASSIGNMENT AND ORDERING WORK
-----------------------------------
For each real pick: find the nearest point on ANY guide contour (within
MAX_ASSIGNMENT_DISTANCE). This assigns the pick to that guide's fibril.
Then, to ORDER the assigned picks along that fibril correctly (not just
group them), each pick is projected onto its guide contour's own
cumulative arc-length - i.e. "how far along this guide path does this
pick's nearest point sit" - and picks are sorted by that value. This
gives a properly ordered trace using your precise pick positions, even
though the guide itself was only sparsely/roughly clicked.

REQUIRES
--------
numpy, pandas, scipy
"""

import numpy as np
import pandas as pd

# ============================== CONFIG ==================================

GUIDE_POINTS_FILE = "guide_points.txt"
TSV_PATH = "coordinates2/Position018_bin4_2_no_mt_overlap.tsv"
OUTPUT_POINTS_FILE = "Position018_guided_trace_points.txt"

MINIMUM_SCORE = 0.27

MAX_ASSIGNMENT_DISTANCE = 15.0

# ==========================================================================


def load_guides(path):
    cols = ["object", "contour", "x", "y", "z"]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)
    df["guide_id"] = df["object"].astype(str) + "_" + df["contour"].astype(str)
    guides = {}
    for gid, group in df.groupby("guide_id", sort=False):
        points = group[["x", "y", "z"]].to_numpy()
        seg_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
        cumulative = np.concatenate([[0], np.cumsum(seg_lengths)])
        guides[gid] = {"points": points, "cumulative": cumulative}
    return guides


def load_picks(tsv_path, min_score):
    df = pd.read_csv(tsv_path, sep="\t")
    return df[df["score"] >= min_score].reset_index(drop=True)


def project_onto_guide(pick_xyz, guide_points, guide_cumulative):
    """Finds the nearest point ON THE GUIDE PATH (interpolating along
    segments, not just the nearest vertex) to pick_xyz, and returns
    (distance, arc_length_position)."""
    best_dist = None
    best_arclen = None

    for i in range(len(guide_points) - 1):
        a, b = guide_points[i], guide_points[i + 1]
        seg_vec = b - a
        seg_len = np.linalg.norm(seg_vec)
        if seg_len == 0:
            continue
        seg_dir = seg_vec / seg_len
        t = np.clip((pick_xyz - a) @ seg_dir, 0, seg_len)
        closest_point = a + t * seg_dir
        dist = np.linalg.norm(pick_xyz - closest_point)
        arclen = guide_cumulative[i] + t

        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_arclen = arclen

    return best_dist, best_arclen


def check_guide_separation(guides, warn_threshold):
    """Warns if any two DIFFERENT guides come suspiciously close to each
    other anywhere along their length - this catches an accidental
    misclick (a guide point placed too near the wrong fibril) before it
    silently causes picks to be misassigned, rather than trusting every
    guide point to be perfectly placed."""
    gids = list(guides.keys())
    warnings = []
    for i in range(len(gids)):
        for j in range(i + 1, len(gids)):
            pts_i = guides[gids[i]]["points"]
            pts_j = guides[gids[j]]["points"]
            for pi in pts_i:
                dists = np.linalg.norm(pts_j - pi, axis=1)
                min_dist = dists.min()
                if min_dist < warn_threshold:
                    warnings.append((gids[i], gids[j], min_dist))
    return warnings


def main():
    guides = load_guides(GUIDE_POINTS_FILE)
    picks = load_picks(TSV_PATH, MINIMUM_SCORE)
    print(f"{len(guides)} guide contours loaded from {GUIDE_POINTS_FILE}")
    print(f"{len(picks)} picks at or above score {MINIMUM_SCORE}")

    separation_warnings = check_guide_separation(guides, warn_threshold=MAX_ASSIGNMENT_DISTANCE)
    if separation_warnings:
        print()
        print("WARNING: some guide points from DIFFERENT fibrils are suspiciously")
        print("close together - this can happen from an accidental misclick, and")
        print("picks near these spots risk being assigned to the wrong fibril:")
        for gid_a, gid_b, dist in separation_warnings:
            print(f"  guide {gid_a} and guide {gid_b} come within {dist:.2f} voxels "
                  f"of each other")
        print("Check these specific guides in 3dmod before trusting the output near")
        print("these locations.")
        print()

    pick_xyz = picks[["x", "y", "z"]].to_numpy()
    assigned_guide = [None] * len(picks)
    assigned_arclen = [None] * len(picks)
    assigned_dist = [None] * len(picks)

    for i, p in enumerate(pick_xyz):
        best_gid, best_dist, best_arclen = None, None, None
        for gid, guide in guides.items():
            dist, arclen = project_onto_guide(p, guide["points"], guide["cumulative"])
            if dist is None:
                continue
            if best_dist is None or dist < best_dist:
                best_dist, best_gid, best_arclen = dist, gid, arclen

        if best_dist is not None and best_dist <= MAX_ASSIGNMENT_DISTANCE:
            assigned_guide[i] = best_gid
            assigned_arclen[i] = best_arclen
            assigned_dist[i] = best_dist

    picks["guide_id"] = assigned_guide
    picks["arclen"] = assigned_arclen
    picks["assign_dist"] = assigned_dist

    n_assigned = picks["guide_id"].notna().sum()
    print(f"{n_assigned} / {len(picks)} picks assigned to a guide "
          f"(within {MAX_ASSIGNMENT_DISTANCE} voxels)")

    assigned = picks.dropna(subset=["guide_id"]).copy()
    assigned = assigned.sort_values(["guide_id", "arclen"])

    output_lines = []
    contour_num = 0
    for gid, group in assigned.groupby("guide_id", sort=False):
        contour_num += 1
        for _, row in group.iterrows():
            output_lines.append(f"1 {contour_num} {row['x']:.3f} {row['y']:.3f} {row['z']:.3f}")
        print(f"  guide {gid}: {len(group)} picks assigned, ordered along guide length")

    with open(OUTPUT_POINTS_FILE, "w") as f:
        f.write("\n".join(output_lines) + "\n")

    print()
    print(f"Wrote {len(output_lines)} ordered points across {contour_num} "
          f"fibrils to {OUTPUT_POINTS_FILE}")
    print()
    print("Next: feed this straight into spline_trace_and_pick.py, or check")
    print("it in 3dmod first (point2model -open, same as always).")


if __name__ == "__main__":
    main()
