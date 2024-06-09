# scDiagnostics: diagnostic functions to assess the quality of cell type annotations in single-cell RNA-seq data

# Overview

The accurate annotation of cell types is a critical step in single-cell RNA-sequencing (scRNA-seq) analysis. While annotation transfer from a reference dataset offers a convenient and automated approach, it can also introduce biases and errors if not performed carefully.

**`scDiagnostics`** is an R package designed to address this challenge by providing a comprehensive set of diagnostic tools for evaluating the quality of cell type annotations in scRNA-seq data. With **`scDiagnostics`**, researchers can systematically assess the compatibility and accuracy of annotations, ensuring reliable and reproducible results in their scRNA-seq analysis workflow.

# Key Features

-   **Systematic Evaluation**: **`scDiagnostics`** offers a suite of functions to evaluate the quality of cell type annotations in scRNA-seq data. This includes assessing dataset imbalance, detecting incompatibilities between query and reference datasets, and identifying potential annotation errors.

-   **Annotation Accuracy**: The package helps ensure accurate cell type annotations by providing tools to assess annotation ambiguity and cluster heterogeneity. Researchers can identify ambiguous or mixed cell populations and adjust annotations accordingly.

-   **Marker Gene Alignment**: Proper marker gene alignment is crucial for consistent cell type definitions. **`scDiagnostics`** enables the comparison of marker genes between query and reference datasets, ensuring that annotations are based on consistent and reliable markers.

-   **Detection of Dataset Incompatibilities**: Incompatibilities between query and reference datasets can lead to erroneous annotations. **`scDiagnostics`** provides functionality to identify and address these incompatibilities, ensuring that annotations can be reliably transferred.

-   **Cluster Heterogeneity Assessment**: The package includes functions to assess the heterogeneity within annotated clusters. This helps researchers identify clusters that may contain multiple cell types or outlier cells, allowing for more accurate and refined annotations.

# Installation

To install the development version of the package from GitHub use the following command:

``` r
BiocManager::install("ccb-hms/scDiagnostics")
```

NOTE: you will need the [remotes](https://cran.r-project.org/web/packages/remotes/index.html) package to install from GitHub.

To build the package vignettes upon installation use:

``` r
BiocManager::install("ccb-hms/scDiagnostics",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

# Usage

To get an overview of the functionality of the package, refer to the [pkgdown website](https://ccb-hms.github.io/scDiagnostics/index.html) for code examples. The complete documentation of each available function in `scDiagnostics`, which includes implementation details and working examples, is available in the [reference tab](https://ccb-hms.github.io/scDiagnostics/reference/index.html).

**`scDiagnostics`** is designed to be user-friendly and integrates seamlessly into any scRNA-seq analysis workflow. By providing robust diagnostic tools, the package helps ensure the accuracy and reliability of cell type annotations, leading to more meaningful and reproducible results.
