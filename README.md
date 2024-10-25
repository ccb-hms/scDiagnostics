# `scDiagnostics`: diagnostic functions to assess the quality of cell type annotations in single-cell RNA-seq data

[![BioC release version](https://www.bioconductor.org/shields/years-in-bioc/scDiagnostics.svg)](https://www.bioconductor.org/packages/release/bioc/html/scDiagnostics.html) 
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

# Overview

The accurate annotation of cell types is a critical step in single-cell RNA-sequencing (scRNA-seq) analysis. While annotation transfer from a reference dataset offers a convenient and automated approach, it can also introduce biases and errors if not performed carefully.

`scDiagnostics` is an R package designed to address this challenge by providing a comprehensive set of diagnostic tools for evaluating the quality of cell type annotations in scRNA-seq data. With `scDiagnostics`, researchers can systematically assess the compatibility and accuracy of annotations, ensuring reliable and reproducible results in their scRNA-seq analysis workflow.

# Installation

## Installation from Bioconductor (Release)

Users interested in using the stable release version of the `scDiagnostics` package: please follow the installation instructions [**here**](https://bioconductor.org/packages/release/bioc/html/scDiagnostics.html). This is the recommended way of installing the package.

## Installation from GitHub (Development)

To install the development version of the package from GitHub use the following command:

``` r
BiocManager::install("ccb-hms/scDiagnostics")
```

To build the package vignettes upon installation use:

``` r
BiocManager::install("ccb-hms/scDiagnostics",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

# Usage

See the [**scDiagnostics website**](https://ccb-hms.github.io/scDiagnostics/) for an overview of the functionality of the package. The individual documentation of each available function incl. usage details and executable examples is available in the [**reference tab**](https://ccb-hms.github.io/scDiagnostics/reference/index.html).

## Key Features

See the [**Getting Started with scDiagnostics**](https://ccb-hms.github.io/scDiagnostics/articles/scDiagnostics.html) vignette for an overview of the main functionality of the package. As listed below, additional vignettes are available for specific diagnostic tasks.

-   [**Visualization of Cell Type Annotations**](https://ccb-hms.github.io/scDiagnostics/articles/VisualizationTools.html): Illustrates the distributions of cell type annotations of the query and reference dataset, allowing the user to identify potential differences in the cell type composition between datasets. Also provides functionality for assessing the impact of frequently used QC criteria on the cell type annotation confidence, allowing the user to identify systematic relationships between QC metrics and the predicted cell type categories.

-   [**Evaluation of Dataset and Marker Gene Alignment**](https://ccb-hms.github.io/scDiagnostics/articles/DatasetAlignment.html): Provides functionality for assessing dataset alignment through quantitative comparison of query-to-reference projections in reduced dimension space. Additional functionality for assessing marker gene expression across datasets allows the user to identify potential misalignments between reference and query on the level of individual genes.

-   [**Detection of Annotation Anomalies**](https://ccb-hms.github.io/scDiagnostics/articles/AnnotationAnomalies.html): Focuses on identifying inconsistencies or anomalies in cell type annotations between the query and reference datasets through comparison of expert annotations with annotations derived from automated methods. By highlighting discrepancies that could be indicative of potential errors, this feature aids in refining and improving the accuracy and reliability of cell type classifications.
