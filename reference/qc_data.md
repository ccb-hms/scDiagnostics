# Quality Control Single-Cell RNA-Seq Dataset

This dataset contains the processed query dataset from the Bunis
haematopoietic stem and progenitor cell data. It has been preprocessed
to include log-normalized counts, QC metrics, SingleR cell type
predictions, and annotation scores.

## Usage

``` r
qc_data
```

## Format

An object of class `SingleCellExperiment` with 500 rows and 750 columns.

## Source

Bunis DG et al. (2021). Single-Cell Mapping of Progressive
Fetal-to-Adult Transition in Human Naive T Cells Cell Rep. 34(1): 108573

## Details

This dataset underwent the following steps:

- Loads the `hpca` reference dataset using `fetchReference` from the
  `celldex` package.

- Loads the QC dataset (Bunis haematopoietic stem and progenitor cell
  data) from Bunis DG et al. (2021).

- Adds QC metrics to the QC dataset using the function
  `addPerCellQCMetrics` from the `scuttle` package.

- Performs log normalization on the QC dataset using the function
  `logNormCounts` from the `scuttle` package.

- Runs SingleR to predict cell types and assigns predicted labels to the
  QC dataset using the function `SingleR` from the `SingleR` package.

- Assigns annotation scores to the QC dataset.

- Selects specific columns (`total`, `SingleR_annotation`,
  `annotation_scores`) from the cell metadata for downstream analysis.

- Selects highly variable genes (HVGs) using the function `getTopHVGs`
  from the `scran` package on the QC dataset.

## References

Bunis DG et al. (2021). Single-Cell Mapping of Progressive
Fetal-to-Adult Transition in Human Naive T Cells Cell Rep. 34(1): 108573

## See also

Use `data("qc_data")` to load and access the resulting quality control
dataset.

## Examples

``` r
# Load and explore the quality control dataset
data("qc_data")
```
