# Mitochondrial Segmentation & Granule Analysis

Analysis of mitochondrial morphology, cristae, and calcium-phosphate granules
from membrane segmentations in cryo-ET tomograms. This README doubles as a
step-by-step usage guide, reconstructed from the actual walkthrough used
when developing this pipeline.

## Overview

```
Multi-label segmentation (MemBrain-seg or Dragonfly multi-ROI export)
              ↓
mito_volume_fill.py   → fills mito interior, saves overlay for QC
              ↓
(validate overlay in Fiji)
              ↓
mito_analyse_and_plot.py  → measurements + interactive HTML dashboard
              ↓
(optional) mito_3D_render.ijm  → rendered 3D figure in Fiji
```

## Step 1 — Get a multi-label segmentation

Segmentation can come from either:

- **MemBrain-seg** (conda environment) run directly on the tomogram, or
- A **multi-ROI segmentation in Dragonfly**, exported as a single-label TIFF
  stack: right-click the ROI → **Export → Export ROI as image stack** →
  choose TIFF, exporting as a binary/label map (not a rendered image).

Either way you should end up with one TIFF where each voxel's value
corresponds to a class (background / membrane / granule).

## Step 2 — Confirm label values

Label values aren't guaranteed to be the same across exports/segmentation
runs, so check before doing anything else:

```bash
python mito_volume_fill.py --check
```

This prints each label value present in the file and its voxel count.
Update `MEMBRANE_LABEL`, `GRANULE_LABEL`, `BACKGROUND_LABEL` at the top of
`mito_volume_fill.py` to match what you see (defaults are `1`, `2`, `3`
respectively — confirmed values from a prior run, but re-check per sample).

## Step 3 — Fill the mitochondrial volume

```bash
python mito_volume_fill.py
```

Set `SEG_PATH` at the top of the script to your segmentation TIFF first.
This produces:
- `mito_filled.tif` — the filled binary volume
- `mito_overlay.tif` — a coloured QC overlay (cyan = filled volume, blue =
  membrane, red = granules)

## Step 4 — Validate visually in Fiji

**Open `mito_overlay.tif`, not `mito_filled.tif`** — the filled volume alone
displays as a plain white/grey mask, which makes it hard to check anything.
The overlay is the one with the coloured channels.

If it opens as grayscale instead of colour, either:
- `Image → Type → RGB Color`, or
- drag-and-drop the file directly into Fiji instead of using `File → Open`

Cyan should fill the full mitochondrial interior, including any granules
that were near the boundary (orphan granule inclusion — see Step 3 logic in
the script).

## Step 5 — Run measurements and generate the dashboard

```bash
python mito_analyse_and_plot.py
```

Set the paths at the top of the script to point at your `mito_filled.tif`
and original segmentation. This produces:
- `granule_measurements.csv` — per-granule size, depth, nearest-neighbour distance
- `single_mito_dashboard.html` — self-contained interactive dashboard:
  spinnable 3D view (real granule shapes via marching cubes), box-whisker
  and dot-strip plots for all metrics, auto-generated interpretation text,
  and summary metric cards

## Step 6 (optional) — 3D rendered figure in Fiji

```
mito_3D_render.ijm
```

Requires the **3D Viewer** plugin (usually bundled with Fiji; if missing:
**Help → Update → Manage Update Sites** → tick "3D Viewer"). Run via
**Plugins → Macros → Run...** and select the file. Update `FILLED_PATH` and
`SEG_PATH` at the top first.

Renders the mito as a translucent light-blue volume with granules as solid
magenta, in Fiji's 3D Viewer, where you can rotate/zoom freely and use
**File → Take Snapshot** for a publication-quality image.

**Parameters to tune if needed:**
- `GRAN_BLUR` (default 4.0) — increase (try 5–8) for smoother-looking
  granule blobs rather than sharp masses
- `GRAN_THRESH` (default 50) — lower if granules disappear after blurring;
  raise if too much noise appears
- Mito transparency (`setTransparency`, default 0.75) — increase toward 0.9
  for a more ghostly outline, decrease toward 0.5 for more solid

## Parameters (confirmed values)

- `MEMBRANE_LABEL = 1`, `GRANULE_LABEL = 2`, `BACKGROUND_LABEL = 3`
- `ORPHAN_GRANULE_MARGIN = 20` voxels (~20 nm)
- `BACKGROUND_DILATION = 15` voxels
- `VOXEL_SIZE_NM = 0.992`
- `mito_3D_render.ijm`: `GRAN_BLUR = 4.0`, `GRAN_THRESH = 50`

## Files in this folder

- [x] `mito_volume_fill.py`
- [x] `mito_analyse_and_plot.py`
- [x] `mito_3D_render.ijm`
- [ ] `requirements.txt` (numpy, scipy, scikit-image, tifffile, matplotlib)

## Notes

- Both Python scripts have hardcoded input paths at the top (e.g. `SEG_PATH`,
  `FILLED_PATH`) rather than command-line arguments — update these per
  sample/tomogram before running, or refactor to accept args if you'll run
  this across many samples.
- `mito_3D_render.ijm` also has a hardcoded absolute path
  (`/ceph/users/loo89671/...`) pointing to your institution's file server —
  worth genericising before making the repo public if you'd rather not have
  that visible.
- If granules appear as disconnected regions that should really be counted
  as one granule, that's a known issue worth checking in
  `mito_analyse_and_plot.py`'s connected-component logic before trusting
  granule counts for the thesis.
