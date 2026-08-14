# Tangle Burden Quantification (Confocal)

FIJI/ImageJ macro pipeline for quantifying tau tangle burden from confocal
z-stacks of organotypic hippocampal slice cultures (hSyn-P301L/S320F tau-FusionRed),
followed by R statistics and publication figures.

## Pipeline

1. **`TangleBurdenMacro.ijm`** — main processing macro:
   - Imports two-channel z-stacks (FusionRed tau channel + Synaptophysin-EGFP)
   - Generates hippocampal ROI mask from the Synaptophysin channel
     (3D Gaussian blur, threshold, hole-filling, dilation/erosion)
   - Segments tau-positive tangles from the tau channel (3D Gaussian blur,
     thresholding)
   - Exports per-sample measurements: hippocampal volume, tangle count,
     total tau-positive volume, burden, density, per-object volume/sphericity
2. **`PostProcessing_TangleBurdenMacro.ijm`** — post-processing of raw macro
   output.
3. **`SampleComparisonStatsMacro.ijm`** — cross-sample comparison.
4. **`burden_analysis.R`** — publication-ready statistics and figures on the
   tangle burden output. Reads `tangles.csv` (expects a `Burden_percent`
   column), and produces:
   - Summary statistics (n, mean, SD, SEM, median, min, max) →
     `burden_summary_statistics.csv`
   - Violin plot with mean/median/SEM overlay → `tangle_burden_violin.pdf`/`.png`
   - Histogram with mean/median lines → `tangle_burden_histogram.pdf`
   - Shapiro-Wilk normality test (printed to console)
   - Formatted publication table (mean ± SD, median, range, 95% CI) →
     `burden_publication_table.csv`
   - Has a commented-out block for grouped comparison (e.g. WT vs Mutant) if
     a `Group` column is added to the input CSV.
   Packages: tidyverse, ggplot2, ggpubr.

## Imaging parameters

- Leica Stellaris 8 confocal, 25x HC FLUOTAR L VISIR dipping objective, CFS air-interface stage
- Voxel size: 0.577 x 0.577 x 0.4985 µm (x, y, z)
- Tau channel Gaussian blur (pre-segmentation): x=0.8, y=0.8, z=0.5 µm
- Volume filters: min 10 µm³ debris filter; max varies by macro (10000 µm³ in
  `TangleBurdenMacro.ijm`, 200 µm³ in `PostProcessing_TangleBurdenMacro.ijm`
  to exclude large outliers — check these match your intent before running)

## Usage

Run in FIJI: `Plugins > Macros > Run...` and select `TangleBurdenMacro.ijm`,
then `PostProcessing_TangleBurdenMacro.ijm`, then
`SampleComparisonStatsMacro.ijm` for the cross-sample summary table.

Then, in R:

```r
# expects tangles.csv in the working directory
Rscript burden_analysis.R
```

Note: `PostProcessing_TangleBurdenMacro.ijm` has a hardcoded
`hippoVol = 5229750` (µm³) — this is sample-specific and needs updating per
tomogram/sample, not left at this value.

## Files in this folder

- [x] `TangleBurdenMacro.ijm`
- [x] `PostProcessing_TangleBurdenMacro.ijm`
- [x] `SampleComparisonStatsMacro.ijm`
- [x] `burden_analysis.R`
- [ ] `renv.lock` (or a plain list of R package versions)

## Notes / things to check before relying on this for the thesis

- The grouped-comparison block (WT vs Mutant) in `burden_analysis.R` is
  commented out — uncomment and adjust if your final analysis needs
  group-level stats, not just pooled.
- `Burden_percent` values are printed to 4 decimal places in the R output —
  sanity-check that precision is appropriate/consistent with the rest of the
  thesis before finalising figures.
