# MERFISH Colitis Query (Day 9 Colitis) Spatial Dataset

This dataset contains a downsampled query (Day 9, DSS-induced colitis at
peak inflammation) subset of the Cadinu et al. (2024) mouse colitis
MERFISH dataset, for use alongside `merfish_reference_data`.

## Usage

``` r
merfish_query_data
```

## Format

An object of class `SpatialExperiment` with 943 rows and 1300 columns.

## Source

Cadinu CA, et al. (2024), processed and hosted at
<https://doi.org/10.5281/zenodo.18274942>.

## Details

See `merfish_reference_data` for the shared processing steps. Cells are
downsampled to at most 260 per cell type. The `tier2_merged` column
retains the original fine-grained label (e.g. distinguishing
"Fibroblast" from "Inflamed Fibroblast"), which is not passed to
scDiagnostics functions directly but is useful for checking which cells
a diagnostic actually flags.

## References

Cadinu CA, et al. (2024). Charting the cellular biogeography in colitis
reveals fibroblast trajectories and coordinated spatial remodeling. Cell
187(8):2010-2028.

## See also

Use `data("merfish_query_data")` to load and access the resulting query
dataset, and `data("merfish_reference_data")` for comparison with the
reference dataset.

## Examples

``` r
# Load and explore the MERFISH query dataset
data("merfish_query_data")
```
