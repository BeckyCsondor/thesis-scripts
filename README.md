# Thesis Scripts

Analysis and visualisation scripts used in my PhD thesis, covering confocal
tau tangle burden quantification, PI/Hoechst viability assays, cryo-electron
tomography of tau fibril organisation, and mitochondrial segmentation/granule
analysis. Pipelines mix FIJI/Python/MATLAB for analysis with R for statistics
and figure generation.

## Repository structure

Organised by analysis pipeline (not by language), since most pipelines
cross multiple languages:

```
01_tangle_burden_confocal/    # FIJI macros for tau burden quantification -> R stats/figures
02_pi_hoechst_viability/      # FIJI batch macro for PI/Hoechst cell viability quantification
03_fibril_organisation/       # IMOD point picks -> spatial stats -> interactive dashboards
04_mitochondrial_segmentation/  # Membrane segmentation -> volume filling -> granule analysis
```

Each folder has its own README describing the pipeline, inputs/outputs, and
how to run it.

## Requirements

- Python 3.x — see notes in each relevant folder (numpy, scipy, matplotlib,
  scikit-learn, scikit-image, tifffile)
- R — packages: tidyverse, ggplot2, ggpubr (see `01_tangle_burden_confocal/`)
- FIJI/ImageJ for the confocal and mito rendering macros
- IMOD (for `.mod` point picks used in fibril organisation)

## Citation

If you use these scripts, please cite [thesis/paper reference — TBD].

## License

TBD — see LICENSE file.
