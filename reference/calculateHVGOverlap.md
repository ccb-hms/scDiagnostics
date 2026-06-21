# Calculate the Overlap Coefficient for Highly Variable Genes

Calculates the overlap coefficient between the sets of highly variable
genes from a reference dataset and a query dataset.

## Usage

``` r
calculateHVGOverlap(reference_genes, query_genes)
```

## Arguments

- reference_genes:

  A character vector of highly variable genes from the reference
  dataset.

- query_genes:

  A character vector of highly variable genes from the query dataset.

## Value

Overlap coefficient, a value between 0 and 1, where 0 indicates no
overlap and 1 indicates complete overlap of highly variable genes
between datasets.

## Details

The overlap coefficient measures the similarity between two gene sets,
indicating how well-aligned reference and query datasets are in terms of
their highly variable genes. This metric is useful in single-cell
genomics to understand the correspondence between different datasets.

The coefficient is calculated using the formula:

\$\$Coefficient(X, Y) = \frac{\|X \cap Y\|}{min(\|X\|, \|Y\|)}\$\$

where X and Y are the sets of highly variable genes from the reference
and query datasets, respectively, \\\|X \cap Y\|\\ is the number of
genes common to both \\X\\ and \\Y\\, and \\min(\|X\|, \|Y\|)\\ is the
size of the smaller set among \\X\\ and \\Y\\.

## References

Luecken et al. Benchmarking atlas-level data integration in single-cell
genomics. Nature Methods, 19:41-50, 2022.

## Author

Anthony Christidis, <anthony-alexander_christidis@hms.harvard.edu>

## Examples

``` r
# Load data
data("reference_data")
data("query_data")

# Selecting highly variable genes
ref_var <- scran::getTopHVGs(reference_data, n = 500)
#> Warning: 'scran::getTopHVGs' is deprecated.
#> Use 'scrapper::chooseHighlyVariableGenes' instead.
#> See help("Deprecated")
#> Warning: 'fitTrendVar' is deprecated.
#> Use 'scrapper::fitVarianceTrend' instead.
#> See help("Deprecated")
#> Warning: 'combineBlocks' is deprecated.
#> See help("Deprecated")
query_var <- scran::getTopHVGs(query_data, n = 500)
#> Warning: 'scran::getTopHVGs' is deprecated.
#> Use 'scrapper::chooseHighlyVariableGenes' instead.
#> See help("Deprecated")
#> Warning: 'fitTrendVar' is deprecated.
#> Use 'scrapper::fitVarianceTrend' instead.
#> See help("Deprecated")
#> Warning: 'combineBlocks' is deprecated.
#> See help("Deprecated")
overlap_coefficient <- calculateHVGOverlap(reference_genes = ref_var,
                                           query_genes = query_var)
overlap_coefficient
#> [1] 0.93
```
