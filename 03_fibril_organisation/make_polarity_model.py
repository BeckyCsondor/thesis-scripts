#!/usr/bin/env python3
"""
make_polarity_model.py

Takes the high-confidence filament results from polarity_analysis.py and
builds a plain point-list text file that groups filaments into two sets
based on the sign of their mean Y-direction (the dominant, most reliable
component we identified - see the direction vectors in the picks CSV).

This text file can then be converted into an IMOD model with two objects
using point2model - one object per polarity group - so you can colour
them differently in 3dmod for visualization.

USAGE
-----
1. Edit the CONFIG paths below if needed (they default to matching
   polarity_analysis.py's OUTPUT_PREFIX).
2. Run:
       python3 make_polarity_model.py
   This writes a text file (default: polarity_groups_points.txt).
3. Convert it to an IMOD model:
       module load imod
       point2model -open polarity_groups_points.txt polarity_groups.mod
   (-open because these are filament traces, not closed shapes)
4. Open it alongside your tomogram:
       3dmod tomograms/Position022_bin4.mrc polarity_groups.mod
5. In 3dmod: Edit > Object > Color, select Object 1, pick a colour (e.g.
   red), then select Object 2, pick a different colour (e.g. green).
   Object 1 = filaments with negative mean Y direction
   Object 2 = filaments with positive mean Y direction

NOTE: this groups by Y-direction sign only, per your request ("more
generally negative or positive Y"). It does NOT re-derive same/opposite
polarity relationships from angle comparisons the way the pairwise CSV
does - a filament running mostly along Y but tilted could, in principle,
sit close to the 0 boundary. Check FILAMENTS_NEAR_ZERO_Y warning below
if this applies to any of your filaments.
"""

import pandas as pd
import subprocess

# ============================== CONFIG ==================================

FILAMENT_DIRECTIONS_CSV = "polarity_withguides/combined_analysis_filaments.csv"
CONTOUR_PATH = "Position002_guided_trace_points.txt"
OUTPUT_POINTS_FILE = "Position002_polarity_groups_points.txt"
OUTPUT_MODEL_FILE = "Position002_polarity_groups.mod"

# REQUIRED when FILAMENT_DIRECTIONS_CSV contains multiple tomograms
# combined (e.g. polarity_and_bundle_analysis.py's output): restricts
# processing to only this tomogram's own filaments, matched by name,
# BEFORE any contour lookup happens. Without this, a filament's local
# contour number (e.g. "1_79") can coincidentally collide with a
# DIFFERENT tomogram's own contour "1_79" - since contour numbering
# restarts from 1 in every tomogram's guide file - silently applying
# the wrong tomogram's polarity classification to this tomogram's real
# points. Set to None only if FILAMENT_DIRECTIONS_CSV already contains
# just one tomogram's filaments.
TOMOGRAM_NAME = "Position002"

# Same confidence thresholds as polarity_analysis.py - only filaments
# meeting both are included in the output model.
MIN_PICKS_PER_FILAMENT = 5
MIN_AGREEMENT = 0.75

# Flag filaments whose mean Y direction is close to zero (i.e. the
# filament isn't running mostly along Y, so a positive/negative split
# may not be meaningful for it).
Y_NEAR_ZERO_THRESHOLD = 0.3

# ==========================================================================


def load_filament_directions(path, min_picks, min_agreement):
    df = pd.read_csv(path)
    confident = df[
        (df["n_picks"] >= min_picks) & (df["internal_agreement"] >= min_agreement)
    ].copy()
    return confident


def load_contours(path):
    cols = ["object", "contour", "x", "y", "z"]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)
    df["filament_id"] = df["object"].astype(str) + "_" + df["contour"].astype(str)
    return df


def main():
    filaments = load_filament_directions(
        FILAMENT_DIRECTIONS_CSV, MIN_PICKS_PER_FILAMENT, MIN_AGREEMENT
    )

    if TOMOGRAM_NAME is not None:
        n_before = len(filaments)
        filaments = filaments[filaments["filament_id"].str.startswith(f"{TOMOGRAM_NAME}::")].copy()
        print(f"Restricted to tomogram '{TOMOGRAM_NAME}': {len(filaments)} / {n_before} "
              f"high-confidence filaments belong to this tomogram")

    contours = load_contours(CONTOUR_PATH)

    print(f"{len(filaments)} high-confidence filaments loaded "
          f"(n_picks >= {MIN_PICKS_PER_FILAMENT}, "
          f"agreement >= {MIN_AGREEMENT})")

    near_zero = filaments[filaments["mean_dir_y"].abs() < Y_NEAR_ZERO_THRESHOLD]
    if not near_zero.empty:
        print()
        print(f"WARNING: {len(near_zero)} filament(s) have a mean Y direction "
              f"close to zero (|Y| < {Y_NEAR_ZERO_THRESHOLD}), meaning they "
              f"don't run mostly along Y. Their group assignment (below) may "
              f"not be meaningful:")
        print(near_zero[["filament_id", "mean_dir_y"]].to_string(index=False))
        print()

    lines = []
    neg_contour_num = 0
    pos_contour_num = 0
    n_neg, n_pos = 0, 0

    for _, row in filaments.iterrows():
        fid = row["filament_id"]
        # filament_id from polarity_analysis.py's multi-tomogram output is
        # prefixed as "TomogramName::object_contour" (e.g.
        # "Position018::1_109"), but the raw contour trace file only has
        # plain "object_contour" IDs (e.g. "1_109") - strip the prefix
        # before looking it up, or the match always fails silently.
        lookup_id = fid.split("::")[-1] if "::" in fid else fid
        points = contours[contours["filament_id"] == lookup_id]
        if points.empty:
            print(f"WARNING: no traced points found for filament {fid} "
                  f"(looked up as '{lookup_id}') in "
                  f"{CONTOUR_PATH} - skipping")
            continue

        if row["mean_dir_y"] < 0:
            obj_num = 1
            neg_contour_num += 1
            contour_num = neg_contour_num
            n_neg += 1
        else:
            obj_num = 2
            pos_contour_num += 1
            contour_num = pos_contour_num
            n_pos += 1

        for _, p in points.iterrows():
            lines.append(f"{obj_num} {contour_num} {p['x']:.3f} {p['y']:.3f} {p['z']:.3f}")

    with open(OUTPUT_POINTS_FILE, "w") as f:
        f.write("\n".join(lines) + "\n")

    print()
    print(f"Object 1 (negative Y direction): {n_neg} filaments")
    print(f"Object 2 (positive Y direction): {n_pos} filaments")
    print(f"Wrote {len(lines)} points to {OUTPUT_POINTS_FILE}")

    print()
    print(f"Running point2model to build {OUTPUT_MODEL_FILE} ...")
    try:
        result = subprocess.run(
            ["point2model", "-open", OUTPUT_POINTS_FILE, OUTPUT_MODEL_FILE],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"Successfully created {OUTPUT_MODEL_FILE}")
        else:
            print(f"point2model failed (exit code {result.returncode}):")
            print(result.stderr)
    except FileNotFoundError:
        print("ERROR: 'point2model' command not found - have you run 'module load imod' "
              "in this terminal session before running this script? point2model is an "
              "IMOD tool, not a Python package, and needs that module loaded to be on PATH.")

    print()
    print("Next: open alongside your tomogram in 3dmod and set colours per object")
    print("via Edit > Object > Color (Object 1, then Object 2).")
    print(f"  module load 3dmod")
    print(f"  3dmod <your_tomogram.mrc> {OUTPUT_MODEL_FILE}")


if __name__ == "__main__":
    main()
