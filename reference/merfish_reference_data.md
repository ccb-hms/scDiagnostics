# MERFISH Colitis Reference (Healthy) Spatial Dataset

This dataset contains a downsampled reference (Day 0, healthy colon)
subset of the Cadinu et al. (2024) mouse colitis MERFISH dataset. It is
a
[SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
on a 943-gene targeted panel, preprocessed to include log-normalized
counts, a cell-type column, and PCA results, restricted to 5 shared cell
types (including fibroblasts).

## Usage

``` r
merfish_reference_data
```

## Format

An object of class `SpatialExperiment` with 943 rows and 1100 columns.

## Source

Cadinu CA, et al. (2024), processed and hosted at
<https://doi.org/10.5281/zenodo.18274942>.

## Details

This dataset underwent the following steps:

- Downloads the processed reference (Day 0) object from Zenodo
  (<https://doi.org/10.5281/zenodo.18274942>), itself derived from the
  Cadinu et al. (2024) MERFISH mouse colitis dataset (`MerfishData`
  package).

- Collapses inflammation-associated variants of a cell state (e.g.
  "Inflamed Fibroblast", present only in the query) into their parent
  lineage for a shared `cell_type_merged` column, while retaining the
  original fine-grained `tier2_merged` label.

- Restricts to 5 shared cell types (Fibroblast, Smooth Muscle,
  Epithelial, Other Immune, Endothelial).

- Downsamples to at most 220 cells per cell type.

- Performs Principal Component Analysis (PCA) using the function
  `runPCA` from the `scater` package on the full 943-gene panel.

## References

Cadinu CA, et al. (2024). Charting the cellular biogeography in colitis
reveals fibroblast trajectories and coordinated spatial remodeling. Cell
187(8):2010-2028.

## See also

Use `data("merfish_reference_data")` to load and access the resulting
reference dataset, and `data("merfish_query_data")` for the
corresponding query (Day 9 colitis) dataset.

## Examples

``` r
# Load and explore the MERFISH reference dataset
data("merfish_reference_data")
```
