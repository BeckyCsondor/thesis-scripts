# Thesis Scripts

Analysis and visualisation scripts used in my PhD thesis, covering confocal
tau tangle burden quantification, PI/Hoechst viability assays, cryo-electron
tomography of tau fibril organisation, and mitochondrial segmentation and
granule analysis. Pipelines mix FIJI and Python for analysis with R for
statistics and figure generation.

## Repository structure

Organised by analysis pipeline rather than by language, since most pipelines
cross several languages:

```
01_tangle_burden_confocal/     FIJI macros for tau burden quantification, R stats and figures
02_pi_hoechst_viability/       FIJI batch macro for PI/Hoechst viability quantification
03_fibril_organisation/        IMOD point picks, spatial statistics, polarity, interactive dashboards
04_mitochondrial_segmentation/ Membrane segmentation, volume filling, granule analysis
05_statistics_and_figures/     R scripts producing the final thesis figures
```

Each folder has its own README describing the pipeline, its inputs and
outputs, and how to run it.

## Requirements

- Python 3.11 — numpy, scipy, pandas, matplotlib, scikit-learn,
  scikit-image, tifffile
- R — tidyverse, ggplot2, ggpubr, dplyr, readr, tidyr
- FIJI/ImageJ for the confocal, viability and mitochondrial rendering macros
- IMOD for the `.mod` point picks used in fibril organisation

## Use of AI assistance

Claude (Anthropic) was used during the development of this work. Its
contribution differed by pipeline:

| Pipeline | Use |
|---|---|
| 01_tangle_burden_confocal | Debugging and refinement of macros written by the author |
| 02_pi_hoechst_viability | Debugging and refinement of macros written by the author |
| 03_fibril_organisation | Scripts substantially written with Claude's assistance, from analysis designs specified by the author |
| 04_mitochondrial_segmentation | Scripts substantially written with Claude's assistance |
| 05_statistics_and_figures | Figure scripts and statistical testing written with Claude's assistance |

In all cases the analyses, parameters and interpretation are the author's
own, and all outputs were checked against the underlying data.

## Citation

Csondor, R. (2026) PhD thesis, University College London.

## License

See LICENSE.
