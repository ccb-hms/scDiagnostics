# Zeisel Brain Reference Single-Cell RNA-Seq Dataset

This dataset contains the processed reference dataset from the Zeisel
mouse brain data. It has been preprocessed to include log-normalized
counts, a cell type column, and PCA results, and is used to benchmark
`detectAnomaly` and `calculateReconstructionError` against known ground
truth.

## Usage

``` r
zeisel_reference_data
```

## Format

An object of class `SingleCellExperiment` with 230 rows and 2103
columns.

## Source

The Zeisel mouse brain dataset, available through the scRNAseq package.

## Details

This dataset underwent the following steps:

- Loads the Zeisel mouse brain dataset from the `scRNAseq` package.

- Performs log normalization using the function `logNormCounts` from the
  `scuttle` package.

- Renames the `level1class` column to `true_cell_type`.

- Divides the data into a 70

- Selects highly variable genes (HVGs) using the function `getTopHVGs`
  from the `scran` package, intersected between the reference and query
  datasets.

- Performs Principal Component Analysis (PCA) on the reference dataset
  using the function `runPCA` from the `scater` package.

## References

Zeisel A, et al. (2015). Cell types in the mouse cortex and hippocampus
revealed by single-cell RNA-seq. Science 347(6226):1138-42.

## See also

Use `data("zeisel_reference_data")` to load and access the resulting
reference dataset, and `data("zeisel_query_data")` for the corresponding
query dataset.

## Examples

``` r
# Load and explore the Zeisel reference dataset
data("zeisel_reference_data")
```
