# Statistics and Figures

R scripts producing the final fibril organisation figures in the thesis.
Each reads a CSV exported from the Python analyses in
`03_fibril_organisation/` and produces a publication-format plot.

All scripts share a common style: Helvetica Neue with a Helvetica fallback,
`ragg` used for rendering where available, and a consistent colour scheme
across related figures.

## Fibril spacing

**`fibril_spacing_overall.R`** — box-and-whisker of overall
nearest-neighbour spacing. Line-to-line distance minus 6 nm fibril width.

**`fibril_spacing_pooled_violin.R`** — violin of the same pooled
distribution, final published version.

**`fibril_spacing_comparison.R`** — pooled nearest neighbour against
orientation-group-restricted nearest neighbour.

**`fibril_spacing_combined.R`** — three-violin figure comparing minimum
approach distance with the sustained parallel gap.

Exclusions applied consistently across these: fibril pairs visually
confirmed in IMOD as fragmented traces of a single fibril rather than
genuine close contacts, upper-IQR outliers, and Position022 throughout
(tracing not considered reliable).

Input: `spacing_overall.csv` (columns: tomogram, gap_nm)

## Alignment

**`pooled_alignment.R`** — S2 nematic order parameter and parallelism index,
both computed within each fibril's own orientation group rather than against
a single whole-tomogram director.

Input: `alignment_pooled.csv` (columns: tomogram, fibril_idx, S2, PI)

## Near-hexagonal coordination centres

**`pooled_nn_clustering.R`** — tests whether annotated nucleation points sit
closer together than equally sized random subsets of fibrils drawn from the
same local density. Real distances are measured only against other centres
within the same view, and the null draws are constructed the same way per
view before pooling, so coordinates from different views are never mixed.

Input: `pooled_nn_data.csv` (columns: group, value)

**`chain_enrichment.R`** — chain membership compared across three groups:
all fibrils as baseline, automated hexagonal candidates selected by a
chain-blind geometric criterion (angular gap standard deviation below 20
degrees), and manually annotated nucleation points.

Input: `validation_summary.csv` (columns: group, n, pct_in_chain,
pvalue_label)

## Feature prevalence statistics

**`prevalence_statistics.py`** — tests the two annotation prevalence tables
reported in the mitochondrial chapter. Feature prevalence is compared
between sample types, and between proximity tiers within the 2xTau sample,
by chi-squared test of independence on the underlying tomogram counts, with
Fisher's exact test for pairwise comparisons and wherever expected cell
counts fall below five.

A separate comparison treats granule-containing mitochondria as a proportion
of all mitochondria observed rather than of all tomograms, which controls for
the differing overall prevalence of mitochondria between groups. Wilson score
intervals are reported for these proportions, since several groups are small.

The counts are held in the script itself rather than read from a file, so it
runs standalone and reproduces the reported values directly.

## Note on test statistics

The underlying statistical tests for these figures were computed in Python
using SciPy. The R scripts read the resulting values and render the figures,
with significance shown as brackets and asterisks rather than raw p-values.
