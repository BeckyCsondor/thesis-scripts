# PI/Hoechst Viability & Tau-Proximity Analysis

Quantifies cell viability in mouse hippocampal organotypic slice cultures
using a propidium iodide (PI) / Hoechst 33342 dual-staining assay, plus a
related analysis of PI+ (dead/dying) nuclei proximity to Tau+ puncta.

## Assay protocol

- Slices washed 3x in PBS
- Staining solution in PBS: PI at 1:100, Hoechst 33342 at 1:1000
- 1 ml applied above the membrane insert, 1 ml below (apical + basolateral)
- Incubated 37°C, 5% CO2, 30 min, protected from light
- Washed 3x in PBS before imaging
- Imaged on Leica Stellaris 8 confocal: AOTF laser transmission 2.5% for
  both PI and Hoechst channels, 10% for tau-FusionRed channel; z-stacks
  tiled in a spiral pattern and stitched for a full hippocampal overview
- Pixel size: 0.433 µm/px

## Viability quantification

**`PI_Hoechst_Viability_Batch.ijm`** — FIJI batch macro. Channel assignment
differs by condition (confirm this matches your acquisition before running):

| Condition | Ch1 | Ch2 | Ch3 |
|---|---|---|---|
| Control | Hoechst | PI (yellow) | — |
| wtTau | Hoechst | Tau (red) | PI (yellow) |
| 2xTau | Tau (red) | Hoechst | PI (yellow) |

Processing:
- Background-subtracts (rolling ball, 50 px)
- Segments Hoechst-positive nuclei (Otsu dark thresholding), size filter
  50–500 px², circularity 0.20–1.00
- Segments PI-positive (dead/dying) nuclei using a manually validated
  threshold (15–255), size filter 10–500 px² (PI+ nuclei are smaller/pyknotic,
  averaging ~34 px²)
- Reports both a segmentation-based and an intensity-based % dead per sample
- Batch-processes multiple open images, exporting a summary CSV

**`viability_plots.R`** — takes the segmentation-based % dead values and
produces:
- Figure A: dot plot, individual slice values + condition mean line
- Figure B: bar graph of mean % dead per condition
- Combined 2-panel figure (via patchwork)
- Saves all three as PDF + PNG

**Note:** the data in `viability_plots.R` is currently hardcoded (pasted-in
values from a prior run: Control 1.6/6.8%, wtTau 7.3/13.0%, 2xTau
36.1/33.5%) rather than reading from the CSV the macro produces — replace
with your actual final `Seg_Pct_Dead` values (or refactor to read the CSV
directly) before treating this as final.

## Tau-proximity analysis

**`NearestNeighbour_TauPI.ijm`** — for a given image, segments Tau+ puncta
(threshold 200–255) and PI+ nuclei (background subtract, Gaussian blur
sigma=1, threshold 15–255, size 10–500 px², circularity 0.20–1.00), then
computes the distance from every PI+ nucleus to its nearest Tau+ punctum.
Prompts for the tau/PI channel window titles at runtime. Prints results to
the FIJI log as CSV-formatted text (`PI_nucleus,Nearest_Tau_dist_um`) — copy
that output out manually into the R script below.

**`nearest_neighbour_plot.R`** — takes those distances (currently hardcoded
per-sample, `s1_distances`/`s2_distances` for 2xTau S1/S2) and produces a
histogram binned by distance, coloured into three proximity zones:
- Same cell (0–5 µm)
- Adjacent cell (5–15 µm) — ~1 cell diameter
- Not associated (>15 µm)

Reports counts/percentages per zone and the median distance as an annotation
on the plot. Saves `nearest_neighbour_TauPI.pdf`/`.png`.

## Confirmed final parameters

| Parameter | Value | Rationale |
|---|---|---|
| PI threshold | 15–255 | Manually validated against real nuclei, not background |
| PI size filter | 10–500 px² | Matches measured PI+ nuclei average (~34 px²) |
| Hoechst size filter | 50–500 px² | Appropriate for larger healthy nuclei |
| Hoechst threshold | Otsu dark | Standard for dense nuclear staining |
| Background subtraction | Rolling ball, 50 px | ~21 µm, larger than any nucleus |
| Circularity filter | 0.20–1.00 | Excludes non-nuclear debris/artifacts |
| Cell diameter (NND analysis) | ~15 µm | Used to define "adjacent cell" zone |

## Usage

```
Fiji: Plugins > Macros > Run... > PI_Hoechst_Viability_Batch.ijm
```
```r
Rscript viability_plots.R
```

For the proximity analysis, per image:
```
Fiji: Plugins > Macros > Run... > NearestNeighbour_TauPI.ijm
# copy the printed distances from the Log window into nearest_neighbour_plot.R
```
```r
Rscript nearest_neighbour_plot.R
```

## Files in this folder

- [x] `PI_Hoechst_Viability_Batch.ijm` — **note:** the file as originally
  exported had `PI_THRESHOLD_LOW = 8`; this has been corrected to `15` here
  to match the confirmed final validated value. Double-check this against
  your own records before trusting it for the thesis.
- [x] `viability_plots.R`
- [x] `NearestNeighbour_TauPI.ijm`
- [x] `nearest_neighbour_plot.R`

## Notes / things to check before relying on this for the thesis

- Both R scripts use hardcoded, pasted-in data rather than reading from a
  CSV — fine for a one-off figure, but a risk if you rerun the macros and
  forget to update the numbers in the R script. Consider refactoring to
  read from CSV directly if you'll regenerate these figures more than once.
- The channel assignment table above differs per condition — this is easy
  to get wrong when batch-processing; worth double-checking image titles
  match the expected condition before running.
- `NearestNeighbour_TauPI.ijm` requires manual copy-paste of results from
  the Fiji Log into the R script — no direct file export. Consider adding a
  `File.saveString(...)` call to the macro if you'll be running this on many
  images, to avoid transcription errors.
