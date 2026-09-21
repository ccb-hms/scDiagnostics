# Zeisel Brain Anomaly Detection Benchmark Results

A precomputed benchmark of `detectAnomaly` (Isolation Forest) and
`calculateReconstructionError` on the full Zeisel mouse brain dataset,
evaluating detection accuracy against known ground truth (a withheld
cell type) across label-noise, class-imbalance, and batch-effect
gradients, plus a hyperparameter grid search. Used in the
`ZeiselBenchmarking` vignette to illustrate how these functions perform
beyond a single worked example.

## Usage

``` r
zeisel_benchmark_results
```

## Format

An object of class `list` of length 3.

## Source

Computed from the Zeisel mouse brain dataset (scRNAseq package) using
[`scDiagnostics::detectAnomaly`](https://ccb-hms.github.io/scDiagnostics/reference/detectAnomaly.md)
and
[`scDiagnostics::calculateReconstructionError`](https://ccb-hms.github.io/scDiagnostics/reference/calculateReconstructionError.md).

## Details

A named list with three data frames:

- `gradients`: AUROC, AUPRC, sensitivity, and specificity for both
  methods across 3 baseline scenarios (a distinct, a related, and a rare
  cell type withheld from the reference) and, for the related and rare
  scenarios, across gradients of label noise (0-40 reference class
  imbalance (300 to 10 cells), and query batch effects (a mean shift
  applied to 20

- `if_tuning`: sensitivity/specificity of `detectAnomaly` across a grid
  of PC subsets, HVG counts, and anomaly thresholds.

- `re_tuning`: sensitivity/specificity of `calculateReconstructionError`
  across a grid of HVG counts, PC subsets, and MAD thresholds.

Computed once, offline, on the full (non-downsampled) Zeisel dataset;
see `inst/script/ZeiselBenchmarkResults.R` for the exact procedure.

## References

Zeisel A, et al. (2015). Cell types in the mouse cortex and hippocampus
revealed by single-cell RNA-seq. Science 347(6226):1138-42.

## See also

Use `data("zeisel_benchmark_results")` to load and access the benchmark
results.

## Examples

``` r
# Load and explore the Zeisel benchmark results
data("zeisel_benchmark_results")
```
