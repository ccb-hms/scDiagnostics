# `scDiagnostics`: diagnostic functions to assess the quality of cell type annotations in single-cell RNA-seq data

# Overview

The accurate annotation of cell types is a critical step in single-cell RNA-sequencing (scRNA-seq) analysis. While annotation transfer from a reference dataset offers a convenient and automated approach, it can also introduce biases and errors if not performed carefully.

`scDiagnostics` is an R package designed to address this challenge by providing a comprehensive set of diagnostic tools for evaluating the quality of cell type annotations in scRNA-seq data. With `scDiagnostics`, researchers can systematically assess the compatibility and accuracy of annotations, ensuring reliable and reproducible results in their scRNA-seq analysis workflow.

# Installation

## Release Version

To install the release version of the package from Bioconductor use the following command:

``` r
BiocManager::install("scDiagnostics")
```

NOTE: you will need the [**BiocManager**](https://cran.r-project.org/web/packages/BiocManager/index.html) package to install from Bioconductor.

To build the package vignettes upon installation use:

``` r
BiocManager::install("scDiagnostics",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

## Development Version

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

You may browse the [**scDiagnostics website**](https://ccb-hms.github.io/scDiagnostics/) website for an overview of the functionality of the package. The individual documentation of each available function in `scDiagnostics`, which includes usage details and executable examples, is available in the [**reference tab**](https://ccb-hms.github.io/scDiagnostics/reference/index.html).

## Key Features

To get an overview of the main functionality of the scDiagnostics package, refer to the [**Getting Started with scDiagnostics**](https://ccb-hms.github.io/scDiagnostics/articles/scDiagnostics.html) vignette. The links below direct you to vignettes that explore functions organized around a common theme or purpose, highlighting how they interact or complement each other in specific contexts.

-   [**Visualization of Cell Type Annotations**](https://ccb-hms.github.io/scDiagnostics/articles/VisualizationTools.html): Illustrates the distributions of cell type annotations of the query and reference dataset, allowing the user to identify potential differences in the cell type composition between datasets.

-   [**Evaluation of QC and Annotation Scores**](https://ccb-hms.github.io/scDiagnostics/articles/QCandAnnotationScores.html): Provides functionality for assessing the impact of frequently used QC criteria on the cell type annotation confidence, allowing the user to identify systematic relationships between QC metrics and the predicted cell type categories.

-   [**Evaluation of Dataset and Marker Gene Alignment**](https://ccb-hms.github.io/scDiagnostics/articles/DatasetAlignment.html): Provides functionality for assessing dataset alignment through quantitative comparison of query-to-reference projections in reduced dimension space. Additional functionality for assessing marker gene expression across datasets allows the user to identify potential misalignments between reference and query on the level of individual genes.

-   [**Statistical Measures to Assess Dataset Alignment**](https://ccb-hms.github.io/scDiagnostics/articles/StatisticalMeasures.html): Utilizes statistical tests and metrics to quantitatively assess the alignment between the query and reference datasets. Key measures include p-values from Hotelling's T-squared test and Cramer's V statistic, which evaluate the degree of similarity or dissimilarity in cell type distributions and principal component projections. These statistical assessments help in determining the robustness of dataset alignment and highlight any significant differences that may impact the reliability of cell type annotations.

-   [**Detection of Annotation Anomalies**](https://ccb-hms.github.io/scDiagnostics/articles/AnnotationAnomalies.html): Focuses on identifying inconsistencies or anomalies in cell type annotations between the query and reference datasets through comparison of expert annotations with annotations derived from automated methods. By highlighting discrepancies that could be indicative of potential errors, this feature aids in refining and improving the accuracy and reliability of cell type classifications.

-   [**Analysis of Distances Between Specific Cells and Cell Populations**](https://ccb-hms.github.io/scDiagnostics/articles/CellDistancesDiagnostics.html): Calculates distances or similarities between individual cells and predefined cell types in both the query and reference datasets. This analysis helps determine how closely each cell in the query dataset matches the cell types defined in the reference dataset. By providing insights into cell type classification accuracy and identifying potential mismatches, this functionality supports more precise annotation and aids in detecting areas that may require further investigation or refinement.
