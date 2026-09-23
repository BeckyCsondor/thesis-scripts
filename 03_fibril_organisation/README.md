# Fibril Organisation Analysis

Analysis of tau fibril spatial organisation within neurofibrillary tangles,
from cryo-ET tomograms. Covers fibril tracing, higher-order architecture,
filament polarity, and the relationship between fibril geometry and
surrounding ultrastructure.

## Tracing pipeline

1. **IMOD** — fibril positions picked as point models (`.mod`) in 3dmod,
   interpolated with `addModPts` to densify points along each contour.
2. **`model2point`** — converts the `.mod` file to a plain text point list
   (`object contour x y z`, space-separated).
3. **`assign_picks_to_manual_guides.py`** — where automated tracing failed
   within dense bundles, a sparse guide contour was traced manually per
   fibril and template-matching picks were snapped onto the nearest guide
   and ordered along its arc length, producing an ordered trace built from
   the real pick positions.
4. **`spline_trace_and_pick.py`** — fits a smooth spline through each
   filament's ordered points and resamples evenly along its arc length,
   producing both a dense curve for visualisation and evenly spaced
   particle positions at a defined number of helical rises. Also writes a
   RELION-style star file with coordinates, helical tube IDs and tangent-
   derived angle priors. The angle prior convention is a geometric
   best effort and should be cross-checked against the RELION-tomo
   conventions in use before being relied on for extraction; the script
   header notes this.

## Architecture analysis

**`fibril_analysis.py`** — the main spatial analysis script. Reads a point
list and computes orientation group clustering (silhouette-optimised k,
with isolated fibrils flagged rather than force-assigned), the S2 nematic
order parameter and parallelism index via the Q-tensor, nearest-neighbour
distances, end-on projections per orientation group, packing-plane
analysis, and chain detection. Outputs a set of summary figures and
`fibril_data.json`.

The script also contains exploratory analyses that were not carried through
to the thesis: psi6 hexagonal packing order at multiple radii, the radial
distribution function, and Voronoi coordination. These remain in the code
but no thesis figure or statistic derives from them.

**`fibril_dashboard.html`** and **`fibril_comparison.html`** — self-contained
interactive dashboards. Open directly in a browser and load one or two
`fibril_data.json` files to explore clusters and statistics visually.

## Polarity

**`polarity_and_bundle_analysis.py`** — assigns template-matching picks to
traced filaments, converts Euler angles to direction vectors, and computes
each filament's mean direction and internal agreement. Pairwise polarity
relationships and spatial bundles are computed within each tomogram
separately and then combined for overall tallies.

**`make_polarity_model.py`** — builds a two-object IMOD model from the
high-confidence filament directions so that the two polarity groups can be
coloured separately in 3dmod for visualisation.

## Fibril geometry relative to surrounding ultrastructure

**`fibril_perturbation_analysis.py`** — for every point along every traced
fibril, computes local curvature, local sustained inter-fibril gap, and
distance to the nearest instance of each annotated ultrastructural feature
(ER and mitochondria by nearest boundary point, ribosomes by centre of
mass).

**`bundle_width_analysis.py`** — adds a local bundle-width metric to that
output. Where the sustained gap measures the distance to the single nearest
qualifying neighbour, bundle width counts how many parallel fibrils travel
alongside each point: neighbours within 8 degrees of parallel, within 100 nm
centreline separation, and offset longitudinally by no more than 20 nm. A
fibril running alone and a fibril within a tightly packed bundle can share
the same nearest-neighbour gap but differ substantially in bundle width.

**`correlate_perturbation_with_features.py`** — runs the statistics on the
above output. For each feature type, it computes a Spearman correlation
between distance-to-feature and each metric (curvature, local sustained
gap, bundle width), and a Mann-Whitney U comparison between the closest and
farthest quartiles of distance. It reads the output of
`bundle_width_analysis.py` by default, and skips any metric not present in
the input, so it can also be run directly on
`fibril_perturbation_analysis.csv`. Results are reported pooled across tomograms, with a
per-tomogram breakdown printed alongside so that any single tomogram driving
an overall result can be identified.

## Analyses performed without a standalone script

Two analyses reported in the thesis were carried out interactively and were
never saved as reusable scripts:

- Manual annotation of fibril chains and near-hexagonal coordination centres
  on end-on projections, using a browser-based annotation tool built for the
  purpose.
- Local radius of curvature for the four representative fibril segments
  passing close to organelles, measured from annotated tomogram screenshots
  by pixel-based line extraction and circular arc fitting.

The statistics arising from the chain and nucleation-centre annotations
(angular gap regularity, matched local null models, cluster bootstrap and
chain enrichment) are reported in the thesis and plotted by the R scripts in
`05_statistics_and_figures/`.

## Usage

```bash
python assign_picks_to_manual_guides.py     # per tomogram
python spline_trace_and_pick.py             # per tomogram
python fibril_analysis.py model2point_output.txt
python polarity_and_bundle_analysis.py      # all tomograms
python make_polarity_model.py               # per tomogram
python fibril_perturbation_analysis.py
python bundle_width_analysis.py
python correlate_perturbation_with_features.py
```

Configuration paths are set at the top of each script.

## Parameters

- `PIXEL_SIZE = 9.92` Å/pixel (2.48 Å at bin 4)
- Fibril width taken as 6 nm throughout when converting centreline
  distances to surface-to-surface gaps
- Sustained gap: pairs within 8 degrees of parallel, at least 50% length
  overlap
- Bundle width: neighbours within 8 degrees of parallel, 100 voxels
  (approximately 99 nm) maximum centreline separation, 20 nm maximum
  longitudinal offset, local axes estimated over a 20 nm arc-length window
- Polarity confidence: at least 5 picks per filament, internal agreement
  at least 0.6
