# Query Single-Cell RNA-Seq Dataset

This dataset contains the processed query dataset from the HeOrganAtlas
dataset for Marrow tissue. It has been preprocessed to include
log-normalized counts, specific metadata columns, annotations based on
SingleR cell type scoring, and PCA, t-SNE, and UMAP results.

## Usage

``` r
query_data
```

## Format

An object of class `SingleCellExperiment` with 392 rows and 503 columns.

## Source

The HeOrganAtlas dataset, available through the scRNAseq package.

## Details

This dataset underwent the following steps:

- Loads the HeOrganAtlas dataset specifically for Marrow tissue from the
  `scRNAseq` package.

- Divides the loaded dataset into a query dataset used for downstream
  analysis.

- Performs log normalization on the query dataset using the function
  `logNormCounts` from the `scuttle` package.

- Selects specific columns (`percent_mito`, `expert_annotation`) from
  the cell metadata for downstream analysis.

- Selects highly variable genes (HVGs) using the function `getTopHVGs`
  from the `scran` package on the query dataset.

- Computes AUC gene set scores using the function `AUCell_calcAUC` from
  the `AUCell` package based on a CD4 T cell signature containing 12
  known CD4 T cell marker genes (IL7R, CCR7, SELL, LEF1, TCF7, LTB,
  KLF2, IL32, CD2, CD3D, CD3E, CD3G) and adds these scores to the query
  dataset as `gene_set_scores`.

- Intersects the highly variable genes between the query and reference
  datasets to obtain common genes for analysis.

- Performs Principal Component Analysis (PCA) on the query dataset using
  the function `runPCA` from the `scater` package.

- Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on the
  query dataset using the function `runTSNE` from the `scater` package.

- Performs Uniform Manifold Approximation and Projection (UMAP) on the
  query dataset using the function `runUMAP` from the `scater` package.

- Adds SingleR annotations (`SingleR_annotation`) and annotation scores
  (`annotation_scores`) to the query dataset using the function
  `SingleR` from the `SingleR` package.

## References

He, et al. (2020). HeOrganAtlas: a comprehensive human organ atlas based
on single-cell RNA sequencing.

## See also

Use `data("query_data")` to load and access the resulting query dataset
and the `data("reference_data")` for comparison with the reference
dataset.

## Examples

``` r
# Load and explore the query dataset
data("query_data")
```
