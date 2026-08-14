# Fibril Organisation Analysis

Analysis of tau fibril spatial organisation within neurofibrillary tangles,
from cryo-ET tomograms.

## Pipeline

1. **IMOD** — fibril positions picked as point models (`.mod`) in 3dmod,
   interpolated with `addModPts`/`addModPts_mce` to densify points along
   each fibril contour.
2. **`model2point`** — converts the `.mod` file to a plain text point list
   (`contour_id x y z`, space-separated).
3. **`fibril_analysis.py`** — the main analysis script. Reads the point list
   and computes:
   - ψ6 (hexagonal packing order) at multiple radii, globally and per
     orientation cluster
   - Clark-Evans-style spacing regularity, RDF (2D radial distribution function)
   - Nearest-neighbour distance (NND), same-cluster vs different-cluster
   - S2 nematic order parameter / alignment (via Q-tensor) + pairwise angle distribution
   - Discrete orientation group clustering (silhouette-optimised k, with
     isolated/ungrouped fibrils flagged rather than force-assigned)
   - True 3D end-on projections per orientation group (IMOD-style view)
   - Packing-plane analysis, row periodicity, row/chain detection
   - Bundle detection and per-bundle end-on analysis
   - Voronoi coordination analysis (global and per-cluster)
   - Outputs: `fibril_analysis_summary.png`, `per_cluster_analysis.png`,
     `alignment_analysis.png`, `hexagonal_3d_analysis.png`,
     `end_on_view.png`, `packing_plane_view.png`, `row_periodicity.png`,
     `row_chains.png`, `orientation_group_endon.png`, and `fibril_data.json`
4. **`fibril_dashboard.html`** — self-contained interactive dashboard;
   drag-and-drop `fibril_data.json` file(s) to explore clusters and stats
   visually. Supports loading multiple datasets for basic cross-dataset viewing.
5. **`fibril_comparison.html`** — dedicated two-dataset side-by-side
   comparison dashboard (distinct from the multi-load view in
   `fibril_dashboard.html`) — load two `fibril_data.json` files into the two
   named loader slots to compare datasets directly.
6. *(planned)* **R script** — to reformat/plot `fibril_data.json` output in
   a specific publication style, matching the pattern of
   `01_tangle_burden_confocal/burden_analysis.R`. Not yet written.

## Status

Python analysis + dashboards: functionally complete. The R statistics/figure
step (see point 6 above) is the remaining piece for this pipeline.

## Usage

```bash
python fibril_analysis.py model2point_output.txt
```

Then open `fibril_dashboard.html` or `fibril_comparison.html` directly in a
browser (no server needed) and load the resulting `fibril_data.json`.

## Files in this folder

- [x] `fibril_analysis.py` — this is the final version (was named
  `fibrilorganisationscript_alignment_endon.py`; renamed here for clarity).
  An earlier, less complete snapshot also existed under the same
  `fibril_analysis.py` name during development — **do not confuse the two**;
  this repo contains only the final, fuller version.
- [x] `fibril_dashboard.html`
- [x] `fibril_comparison.html`
- [ ] R script for publication-style figures (planned, not yet written)
- [ ] `requirements.txt` (numpy, scipy, matplotlib, scikit-learn)

## Parameters worth double-checking before finalising for the thesis

- `PIXEL_SIZE = 9.92` Å/pixel (2.48 Å × bin4) — confirm this matches every
  tomogram's actual binning; other scripts in this repo use different voxel
  sizes for different datasets.
- `PSI6_RADII = [300, 400, 500, 700]` Å — chosen relative to a specific
  dataset's mean NND (42.6 nm); re-derive if applying to a dataset with
  substantially different fibril spacing.
