# Zeisel Brain Query Single-Cell RNA-Seq Dataset

This dataset contains the processed query dataset from the Zeisel mouse
brain data, for use alongside `zeisel_reference_data`.

## Usage

``` r
zeisel_query_data
```

## Format

An object of class `SingleCellExperiment` with 230 rows and 902 columns.

## Source

The Zeisel mouse brain dataset, available through the scRNAseq package.

## Details

See `zeisel_reference_data` for the processing steps shared by both
datasets.

## References

Zeisel A, et al. (2015). Cell types in the mouse cortex and hippocampus
revealed by single-cell RNA-seq. Science 347(6226):1138-42.

## See also

Use `data("zeisel_query_data")` to load and access the resulting query
dataset, and `data("zeisel_reference_data")` for comparison with the
reference dataset.

## Examples

``` r
# Load and explore the Zeisel query dataset
data("zeisel_query_data")
```
