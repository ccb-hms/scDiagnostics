# COVID-19 PBMC Reference (Healthy) Single-Cell RNA-Seq Dataset

This dataset contains a downsampled reference (healthy donor) subset of
the Stephenson et al. (2021) COVID-19 PBMC atlas. It has been
preprocessed to include log-normalized counts, an author-provided
cell-type column, and PCA results, restricted to 5 shared cell types
(including CD14 monocytes) and a gene panel that always includes the
Yoshida et al. interferon-response signature.

## Usage

``` r
covid_reference_data
```

## Format

An object of class `SingleCellExperiment` with 414 rows and 900 columns.

## Source

Stephenson E, et al. (2021), processed and hosted at
<https://doi.org/10.5281/zenodo.18274942>.

## Details

This dataset underwent the following steps:

- Downloads the processed reference (healthy) object from Zenodo
  (<https://doi.org/10.5281/zenodo.18274942>), itself derived from the
  Stephenson et al. (2021) COVID-19 PBMC atlas (CZI CELLxGENE).

- Restricts to 5 shared cell types (CD14 monocytes, CD4 T, CD8 T, B
  cells, NK_16hi) using the `author_cell_type_merged` column.

- Downsamples to at most 180 cells per cell type.

- Selects highly variable genes using the function `getTopHVGs` from the
  `scran` package, intersected between reference and query, always
  including the 25-gene Yoshida et al. (2019) interferon-response
  signature.

- Performs Principal Component Analysis (PCA) using the function
  `runPCA` from the `scater` package.

## References

Stephenson E, et al. (2021). Single-cell multi-omics analysis of the
immune response in COVID-19. Nature Medicine 27:904-916.

## See also

Use `data("covid_reference_data")` to load and access the resulting
reference dataset, and `data("covid_query_data")` for the corresponding
query (severe COVID-19) dataset.

## Examples

``` r
# Load and explore the COVID-19 reference dataset
data("covid_reference_data")
```
