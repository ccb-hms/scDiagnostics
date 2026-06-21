# Reference Single-Cell RNA-Seq Dataset

This dataset contains the processed reference dataset from the
HeOrganAtlas dataset for Marrow tissue. It has been preprocessed to
include log-normalized counts, specific metadata columns, and PCA,
t-SNE, and UMAP results.

## Usage

``` r
reference_data
```

## Format

An object of class `SingleCellExperiment` with 392 rows and 1500
columns.

## Source

The HeOrganAtlas dataset, available through the scRNAseq package.

## Details

This dataset underwent the following steps:

- Loads the HeOrganAtlas dataset specifically for Marrow tissue from the
  `scRNAseq` package.

- Divides the loaded dataset into a reference dataset used for
  downstream analysis.

- Performs log normalization on the reference dataset using the function
  `logNormCounts` from the `scuttle` package.

- Selects the column `expert_annotation`) from the cell metadata for
  downstream analysis.

- Selects highly variable genes (HVGs) using the function `getTopHVGs`
  from the `scran` package on the reference dataset.

- Performs Principal Component Analysis (PCA) on the reference dataset
  using the function `runPCA` from the `scater` package.

- Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on the
  reference dataset using the function `runTSNE` from the `scater`
  package.

- Performs Uniform Manifold Approximation and Projection (UMAP) on the
  reference dataset using the function `runUMAP` from the `scater`
  package.

## References

He, et al. (2020). HeOrganAtlas: a comprehensive human organ atlas based
on single-cell RNA sequencing.

## See also

Use `data("reference_data")` to load and access the resulting reference
dataset.

## Examples

``` r
# Load and explore the reference dataset
data("reference_data")
```
