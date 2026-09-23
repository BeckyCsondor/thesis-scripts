#!/usr/bin/env python3
"""
spline_trace_and_pick.py

Takes ordered per-filament point traces (e.g. from trace_from_picks.py,
or a model2point export of manual IMOD traces) and:

  1. Fits a smooth spline through each filament's points - giving you a
     clean curve for visualization (e.g. thesis figures showing "visual
     segmentation" of traced fibrils), rather than the raw, slightly
     jagged sequence of individual picks.

  2. Resamples evenly along that spline's ARC LENGTH (not evenly in
     the spline's internal parameter, which is NOT the same thing) at
     a physical spacing you specify - giving you new particle positions
     for extraction, spaced by a defined number of helical rises.

  3. Computes a tangent-direction vector at each resampled position
     (from the spline's derivative) and writes a RELION-style star file
     with coordinates, helical tube IDs, track length, and angle
     priors derived from that tangent.

IMPORTANT - ANGLE PRIOR CONVENTION IS A BEST EFFORT, VERIFY BEFORE USE
------------------------------------------------------------------------
This script computes rlnAngleTiltPrior/rlnAnglePsiPrior as the spherical
polar/azimuthal angle of each tangent vector (tilt = angle from Z axis,
psi = azimuth in the XY plane). This is a defensible general approach
for a fully 3D tomographic tangent (unlike the SPA-only convention of
fixing tilt=90, which assumes the filament lies flat in a 2D image),
but it has NOT been validated against your specific RELION-tomo /
Warp pipeline's exact conventions.

A published in situ cryo-ET filament pipeline using a very similar
Warp/IMOD/RELION workflow to yours built a dedicated, tested tool for
this exact conversion (from imodinfo output to RELION helical priors):
    https://github.com/jjenkins01/model2helicalpriors
It's worth cross-checking (or just using directly) rather than trusting
this script's angle priors blindly for a real extraction run - treat
what's here as a solid geometric starting point, not gospel.

Since rot (rotation about the helical axis itself) can't be determined
from tracing alone, it's intentionally left as a free parameter here -
that's what subtomogram alignment is for. rlnAnglePsiFlipRatio is set
to 0.5 (standard for filaments where each segment's absolute polarity
along the tube isn't independently known) - once you've run your
polarity_analysis.py pipeline on a given filament with a confident
result, you could tighten psi/flip ratio for that filament specifically
rather than leaving it at the generic default.

REQUIRES
--------
scipy, numpy, pandas (already in your polarity_analysis conda env)
"""

import numpy as np
import pandas as pd
from scipy.interpolate import splprep, splev

# ============================== CONFIG ==================================

INPUT_POINTS_FILE = "auto_trace_points.txt"   # object contour x y z, ordered per contour
OUTPUT_PREFIX = "Position022_bin2_spline"

TOMOGRAM_NAME = "Position022_bin2"   # goes into rlnTomoName
PIXEL_SIZE_ANGSTROM = 4.96           # bin2 pixel size

# Desired inter-particle spacing along the fibril, in helical rises.
RISE_ANGSTROM = 4.98
N_RISES_PER_STEP = 10   # -> ~95% box overlap for a ~992 Angstrom box, see chat

# How finely to sample the spline for the VISUALIZATION output (dense
# curve for figures) - in voxels between points, not the same as the
# particle-picking spacing above.
VIS_STEP_ANGSTROM = 5.0

# Spline smoothing factor (scipy's `s` parameter in splprep). 0 = spline
# passes exactly through every input point (can look jagged if picks are
# noisy). Higher values smooth more aggressively. A common starting
# heuristic is roughly the number of points in the trace; adjust and
# re-look at the visualization output if it looks too wiggly or too
# over-smoothed/cut-corners.
SPLINE_SMOOTHING_PER_POINT = 1.0

# Traces with fewer than this many points are skipped (can't fit a
# meaningful spline / not a real filament).
MIN_POINTS_FOR_SPLINE = 4

# ==========================================================================


def load_ordered_traces(path):
    cols = ["object", "contour", "x", "y", "z"]
    df = pd.read_csv(path, sep=r"\s+", header=None, names=cols)
    df["filament_id"] = df["object"].astype(str) + "_" + df["contour"].astype(str)
    traces = {}
    for fid, group in df.groupby("filament_id", sort=False):
        traces[fid] = group[["x", "y", "z"]].to_numpy()
    return traces


def fit_spline(points_xyz, smoothing):
    """Fits a smooth 3D spline through ordered points. Returns
    (tck, u) as used by scipy.interpolate.splev, or None if there
    weren't enough points for a stable fit."""
    n = len(points_xyz)
    if n < MIN_POINTS_FOR_SPLINE:
        return None
    k = min(3, n - 1)  # cubic if possible, lower degree for short traces
    try:
        tck, u = splprep(points_xyz.T, s=smoothing, k=k)
    except Exception as e:
        print(f"  WARNING: spline fit failed ({e}) - skipping this trace")
        return None
    return tck


def arc_length_lookup_table(tck, n_samples=2000):
    """Builds a fine-grained table mapping spline parameter u -> physical
    arc length (in voxels) travelled from the start, by densely sampling
    and summing Euclidean distances. Needed because evenly-spaced-in-u
    is NOT the same as evenly-spaced-in-physical-distance."""
    u_fine = np.linspace(0, 1, n_samples)
    pts = np.array(splev(u_fine, tck)).T
    seg_lengths = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cumulative = np.concatenate([[0], np.cumsum(seg_lengths)])
    return u_fine, cumulative


def resample_by_arc_length(tck, u_fine, cumulative, step_voxels):
    """Returns (points, tangents, track_lengths_voxels) evenly spaced by
    physical arc length along the spline, at the given step."""
    total_length = cumulative[-1]
    if total_length <= 0:
        return np.empty((0, 3)), np.empty((0, 3)), np.empty(0)

    target_lengths = np.arange(0, total_length, step_voxels)
    target_us = np.interp(target_lengths, cumulative, u_fine)

    points = np.array(splev(target_us, tck)).T
    derivatives = np.array(splev(target_us, tck, der=1)).T
    norms = np.linalg.norm(derivatives, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    tangents = derivatives / norms

    return points, tangents, target_lengths


def tangent_to_tilt_psi(tangent_xyz):
    """Spherical polar/azimuthal angle of a 3D tangent vector.
    tilt: angle from +Z axis (0-180 deg). psi: azimuth in XY plane
    (-180 to 180 deg). See the module docstring's caveat about this
    convention needing validation against your specific pipeline."""
    x, y, z = tangent_xyz
    tilt = np.degrees(np.arccos(np.clip(z, -1.0, 1.0)))
    psi = np.degrees(np.arctan2(y, x))
    return tilt, psi


def main():
    traces = load_ordered_traces(INPUT_POINTS_FILE)
    print(f"Loaded {len(traces)} traces from {INPUT_POINTS_FILE}")

    step_voxels = (RISE_ANGSTROM * N_RISES_PER_STEP) / PIXEL_SIZE_ANGSTROM
    vis_step_voxels = VIS_STEP_ANGSTROM / PIXEL_SIZE_ANGSTROM
    print(f"Particle spacing: {N_RISES_PER_STEP} rises x {RISE_ANGSTROM} A "
          f"= {RISE_ANGSTROM * N_RISES_PER_STEP:.2f} A = {step_voxels:.2f} px "
          f"at {PIXEL_SIZE_ANGSTROM} A/px")
    print()

    vis_lines = []
    star_rows = []
    vis_contour_num = 0
    n_traces_used = 0
    n_traces_skipped = 0

    for fid, points in traces.items():
        smoothing = SPLINE_SMOOTHING_PER_POINT * len(points)
        tck = fit_spline(points, smoothing)
        if tck is None:
            n_traces_skipped += 1
            continue
        n_traces_used += 1

        u_fine, cumulative = arc_length_lookup_table(tck)
        total_length_voxels = cumulative[-1]

        # --- dense points for visualization ---
        vis_pts, _, _ = resample_by_arc_length(tck, u_fine, cumulative, vis_step_voxels)
        vis_contour_num += 1
        for p in vis_pts:
            vis_lines.append(f"1 {vis_contour_num} {p[0]:.3f} {p[1]:.3f} {p[2]:.3f}")

        # --- evenly spaced particle picks ---
        pick_pts, pick_tangents, track_lengths_voxels = resample_by_arc_length(
            tck, u_fine, cumulative, step_voxels
        )
        for p, t, track_len_vox in zip(pick_pts, pick_tangents, track_lengths_voxels):
            tilt, psi = tangent_to_tilt_psi(t)
            star_rows.append({
                "rlnTomoName": TOMOGRAM_NAME,
                "rlnCoordinateX": p[0],
                "rlnCoordinateY": p[1],
                "rlnCoordinateZ": p[2],
                "rlnHelicalTubeID": fid,
                "rlnHelicalTrackLengthAngst": track_len_vox * PIXEL_SIZE_ANGSTROM,
                "rlnAngleTiltPrior": tilt,
                "rlnAnglePsiPrior": psi,
                "rlnAnglePsiFlipRatio": 0.5,
            })

        print(f"  {fid}: {len(points)} input picks -> spline length "
              f"{total_length_voxels:.1f} px -> {len(pick_pts)} particles, "
              f"{len(vis_pts)} visualization points")

    with open(f"{OUTPUT_PREFIX}_visualization_points.txt", "w") as f:
        f.write("\n".join(vis_lines) + "\n")

    star_df = pd.DataFrame(star_rows)
    # RELION-tomo requires an integer-like tube ID per tomogram, not a
    # string - remap each filament_id to a small sequential integer.
    tube_id_map = {fid: i + 1 for i, fid in enumerate(star_df["rlnHelicalTubeID"].unique())}
    star_df["rlnHelicalTubeID"] = star_df["rlnHelicalTubeID"].map(tube_id_map)

    star_path = f"{OUTPUT_PREFIX}_particles.star"
    with open(star_path, "w") as f:
        f.write("data_particles\n\nloop_\n")
        for i, col in enumerate(star_df.columns, start=1):
            f.write(f"_{col} #{i}\n")
        for _, row in star_df.iterrows():
            f.write(" ".join(
                f"{v:.6f}" if isinstance(v, float) else str(v) for v in row
            ) + "\n")

    print()
    print(f"{n_traces_used} traces processed, {n_traces_skipped} skipped "
          f"(fewer than {MIN_POINTS_FOR_SPLINE} points)")
    print(f"Saved dense visualization curve points to "
          f"{OUTPUT_PREFIX}_visualization_points.txt")
    print(f"  -> convert with: point2model -open "
          f"{OUTPUT_PREFIX}_visualization_points.txt {OUTPUT_PREFIX}_spline.mod")
    print(f"Saved {len(star_df)} particles to {star_path}")
    print()
    print("Reminder: angle priors in this star file are a best-effort")
    print("geometric convention - verify against model2helicalpriors or")
    print("your RELION-tomo pipeline's expected convention before using")
    print("for a real extraction run (see this script's header).")


if __name__ == "__main__":
    main()
