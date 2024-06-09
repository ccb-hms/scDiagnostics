---
editor_options: 
  markdown: 
    wrap: 72
---

# scDiagnostics: diagnostic functions to assess the quality of cell type annotations in single-cell RNA-seq data

# Purpose

Annotation transfer from a reference dataset for the cell type
annotation of a new query single-cell RNA-sequencing (scRNA-seq)
experiment is an integral component of the typical analysis workflow.
The approach provides a fast, automated, and reproducible alternative to
the manual annotation of cell clusters based on marker gene expression.
However, dataset imbalance and undiagnosed incompatibilities between
query and reference dataset can lead to erroneous annotation and distort
downstream applications.

The `scDiagnostics` package provides functionality for the systematic
evaluation of cell type assignments in scRNA-seq data. `scDiagnostics`
offers a suite of diagnostic functions to assess whether both (query and
reference) datasets are aligned, ensuring that annotations can be
transferred reliably. `scDiagnostics` also provides functionality to
assess annotation ambiguity, cluster heterogeneity, and marker gene
alignment. The implemented functionality helps researchers to determine
how accurately cells from a new scRNA-seq experiment can be assigned to
known cell types.

# Installation

To install the development version of the package from GitHub use the
following command:

``` r
BiocManager::install("ccb-hms/scDiagnostics")
```

NOTE: you will need the
[remotes](https://cran.r-project.org/web/packages/remotes/index.html)
package to install from GitHub.

To build the package vignettes upon installation use:

``` r
BiocManager::install("ccb-hms/scDiagnostics",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

