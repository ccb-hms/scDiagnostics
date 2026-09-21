# COVID-19 PBMC Query (Severe COVID-19) Single-Cell RNA-Seq Dataset

This dataset contains a downsampled query (severe COVID-19 donor) subset
of the Stephenson et al. (2021) COVID-19 PBMC atlas, for use alongside
`covid_reference_data`. Cell types are annotated using Azimuth reference
mapping (`azimuth_celltype_l1_merged`).

## Usage

``` r
covid_query_data
```

## Format

An object of class `SingleCellExperiment` with 414 rows and 1330
columns.

## Source

Stephenson E, et al. (2021), processed and hosted at
<https://doi.org/10.5281/zenodo.18274942>.

## Details

See `covid_reference_data` for the shared processing steps. CD14
monocytes, the focus of the accompanying case study vignette, are
downsampled to at most 450 cells (other cell types to at most 220) so
the interferon-activated subset described in the manuscript stays
well-represented.

## References

Stephenson E, et al. (2021). Single-cell multi-omics analysis of the
immune response in COVID-19. Nature Medicine 27:904-916.

## See also

Use `data("covid_query_data")` to load and access the resulting query
dataset, and `data("covid_reference_data")` for comparison with the
reference dataset.

## Examples

``` r
# Load and explore the COVID-19 query dataset
data("covid_query_data")
```
