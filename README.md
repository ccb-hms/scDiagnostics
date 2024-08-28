# `scDiagnostics`: diagnostic functions to assess the quality of cell type annotations in single-cell RNA-seq data

# Overview

The accurate annotation of cell types is a critical step in single-cell RNA-sequencing (scRNA-seq) analysis. While annotation transfer from a reference dataset offers a convenient and automated approach, it can also introduce biases and errors if not performed carefully.

`scDiagnostics` is an R package designed to address this challenge by providing a comprehensive set of diagnostic tools for evaluating the quality of cell type annotations in scRNA-seq data. With `scDiagnostics`, researchers can systematically assess the compatibility and accuracy of annotations, ensuring reliable and reproducible results in their scRNA-seq analysis workflow.

# Installation

To install the development version of the package from GitHub use the following command:

``` r
devtools::install_github("ccb-hms/scDiagnostics")
```

NOTE: you will need the [remotes](https://cran.r-project.org/web/packages/remotes/index.html) package to install from GitHub.

To build the package vignettes upon installation use:

``` r
devtools::install_github("ccb-hms/scDiagnostics",
                          build_vignettes = TRUE,
                          dependencies = TRUE)
```

# Usage

You may browse the [scDiagnostics website](https://ccb-hms.github.io/scDiagnostics/) website for an overview of the functionality of the package. The complete documentation of each available function in `scDiagnostics`, which includes implementation details and working examples, is available in the [reference tab](https://ccb-hms.github.io/scDiagnostics/reference/index.html).

## Key Features

To get an overview of the main functionality of the scDiagnostics package, refer to the [Getting Started with scDiagnostics](https://ccb-hms.github.io/scDiagnostics/articles/scDiagnostics.html) vignette. The links below direct you to vignettes that explore functions organized around a common theme or purpose, highlighting how they interact or complement each other in specific contexts.

-   [Visualization of Cell Type Annotations](https://ccb-hms.github.io/scDiagnostics/articles/VisualizationTools.html): Provides graphical representations of cell type annotations from both the query and reference datasets. This visualization enables a direct comparison of how cell types are distributed and annotated in each dataset, highlighting any significant differences or similarities in the cell type classification across datasets.

-   [Evaluation of QC and Annotation Scores](https://ccb-hms.github.io/scDiagnostics/articles/QCandAnnotationScores.html): Assesses the quality control (QC) metrics and annotation scores for cells in the query and reference datasets. This involves evaluating various QC parameters and comparing annotation scores to ensure consistency and accuracy. By analyzing these scores, researchers can identify potential issues in data quality and annotation reliability, facilitating improvements in data preprocessing and annotation accuracy. This evaluation helps in ensuring that both datasets meet the necessary quality standards for robust and reliable cell type analysis.

-   [Evaluation of Dataset and Marker Gene Alignment](https://ccb-hms.github.io/scDiagnostics/articles/DatasetAlignment.html): Evaluates the alignment of datasets and marker gene expressions between the query and reference datasets. This often involves projecting query data onto principal components derived from the reference dataset and comparing these projections. In addition, by examining how well the marker genes align across datasets, researchers can identify misalignments or deviations that may indicate batch effects or inconsistencies in cell type annotations, facilitating more accurate dataset integration and interpretation.

-   [Statistical Measures to Assess Dataset Alignment](https://ccb-hms.github.io/scDiagnostics/articles/StatisticalMeasures.html): Utilizes statistical tests and metrics to quantitatively assess the alignment between the query and reference datasets. Key measures include p-values from Hotelling's T-squared test and Cramer's V statistic, which evaluate the degree of similarity or dissimilarity in cell type distributions and principal component projections. These statistical assessments help in determining the robustness of dataset alignment and highlight any significant differences that may impact the reliability of cell type annotations.

-   [Detection of Annotation Anomalies](https://ccb-hms.github.io/scDiagnostics/articles/AnnotationAnomalies.html): Focuses on identifying inconsistencies or anomalies in cell type annotations between the query and reference datasets. This involves comparing expert-generated annotations with those produced by automated methods. By highlighting potential errors or discrepancies, this feature aids in refining and improving the accuracy of cell type classifications, ensuring that annotation results are reliable and accurate.

-   [Analysis of Distances Between Specific Cells and Cell Populations](https://ccb-hms.github.io/scDiagnostics/articles/CellDistancesDiagnostics.html): Calculates distances or similarities between individual cells and predefined cell types in both the query and reference datasets. This analysis helps determine how closely each cell in the query dataset matches the cell types defined in the reference dataset. By providing insights into cell type classification accuracy and identifying potential mismatches, this functionality supports more precise annotation and aids in detecting areas that may require further investigation or refinement.
