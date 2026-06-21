# scDiagnostics: Single-Cell Diagnostics Package

\`scDiagnostics\` is a comprehensive toolkit designed for the analysis
and diagnostics of single-cell RNA sequencing (scRNA-seq) data. This
package provides functionalities for comparing principal components,
visualizing canonical correlation analysis (CCA) outputs, and plotting
cell type-specific MDS and PCA projections.

## Details

The package includes the following key functionalities, organized by
their specific purposes:

## Visualization of Cell Type Annotations

Functions for visualizing differences between query and reference
datasets across multiple cell types.

- [`boxplotPCA`](https://ccb-hms.github.io/scDiagnostics/reference/boxplotPCA.md):
  Boxplots of PCA scores for cell types.

- [`calculateDiscriminantSpace`](https://ccb-hms.github.io/scDiagnostics/reference/calculateDiscriminantSpace.md):
  Calculates discriminant space, with a plot method for visualization.

- [`plotCellTypeMDS`](https://ccb-hms.github.io/scDiagnostics/reference/plotCellTypeMDS.md):
  Creates MDS plots for cell types using query and reference datasets.

- [`plotCellTypePCA`](https://ccb-hms.github.io/scDiagnostics/reference/plotCellTypePCA.md):
  Plots principal components for different cell types.

## Evaluation of Dataset Alignment

Functions for visualizing differences between query and reference
datasets for a specific cell type.

- [`calculateGraphIntegration`](https://ccb-hms.github.io/scDiagnostics/reference/calculateGraphIntegration.md):
  Calculate a graph network and form communities of cells.

- [`calculateWassersteinDistance`](https://ccb-hms.github.io/scDiagnostics/reference/calculateWassersteinDistance.md):
  Wasserstein distance for different cell types.

- [`comparePCA`](https://ccb-hms.github.io/scDiagnostics/reference/comparePCA.md):
  Compares PCA results, with a plot method for visualization.

- [`comparePCASubspace`](https://ccb-hms.github.io/scDiagnostics/reference/comparePCASubspace.md):
  Compares PCA subspace, with a plot method for visualization.

- [`plotPairwiseDistancesDensity`](https://ccb-hms.github.io/scDiagnostics/reference/plotPairwiseDistancesDensity.md):
  Plots the density of pairwise distances.

## Calculation of Statistical Measures to Compare Two Datasets

Functions to compute statistical measures to compare two datasets.

- [`calculateAveragePairwiseCorrelation`](https://ccb-hms.github.io/scDiagnostics/reference/calculateAveragePairwiseCorrelation.md):
  Calculates average pairwise correlation, with a plot method for
  visualization.

- [`calculateCramerPValue`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCramerPValue.md):
  Calculates the p-value using Cramer's test.

- [`calculateHotellingPValue`](https://ccb-hms.github.io/scDiagnostics/reference/calculateHotellingPValue.md):
  Calculates the p-value using Hotelling's T-squared test.

- [`calculateMMDPValue`](https://ccb-hms.github.io/scDiagnostics/reference/calculateMMDPValue.md):
  Calculates the p-value using maximum mean discrepancy.

- [`regressPC`](https://ccb-hms.github.io/scDiagnostics/reference/regressPC.md):
  Performs regression on principal components, with a plot method for
  visualization.

## Evaluation of Marker Gene Alignment

Functions for calculating overlap measures of genes between two
datasets.

- [`calculateHVGOverlap`](https://ccb-hms.github.io/scDiagnostics/reference/calculateHVGOverlap.md):
  Calculates overlap of highly variable genes (HVG) between datasets.

- [`calculateVarImpOverlap`](https://ccb-hms.github.io/scDiagnostics/reference/calculateVarImpOverlap.md):
  Calculates overlap of variable importance measures between datasets.

- [`compareMarkers`](https://ccb-hms.github.io/scDiagnostics/reference/compareMarkers.md):
  Compares marker genes between datasets.

## Visualization of Marker Expressions

Functions for to visualize and compare the expression of markers between
a reference and a query dataset.

- [`plotGeneExpressionDimred`](https://ccb-hms.github.io/scDiagnostics/reference/plotGeneExpressionDimred.md):
  Plots gene expression in a dimensionality reduction space.

- [`plotMarkerExpression`](https://ccb-hms.github.io/scDiagnostics/reference/plotMarkerExpression.md):
  Plots marker expression levels.

## Anomaly Detection (Global and Cell Type-Specific)

Functions for detecting anomalies at both the global and cell
type-specific levels.

- [`detectAnomaly`](https://ccb-hms.github.io/scDiagnostics/reference/detectAnomaly.md):
  Detects anomalies in the data, with a plot method for visualization.

- [`calculateCellSimilarityPCA`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellSimilarityPCA.md):
  Calculates cell similarity in PCA space, with a plot method for
  visualization.

## Calculation of Distances Between Specific Cells and Cell Populations

Functions for calculating distances between specific cells and cell
populations.

- [`calculateCellDistances`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellDistances.md):
  Calculates distances between cells, with a plot method for
  visualization.

- [`calculateCellDistancesSimilarity`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellDistancesSimilarity.md):
  Calculates similarity based on cell distances, with a plot method for
  visualization.

## Visualization of QC and Annotation Scores

Functions for visualizing quality control (QC) metrics or other
characteristics of the data.

- [`histQCvsAnnotation`](https://ccb-hms.github.io/scDiagnostics/reference/histQCvsAnnotation.md):
  Plots histograms of QC metrics versus annotations.

- [`plotGeneSetScores`](https://ccb-hms.github.io/scDiagnostics/reference/plotGeneSetScores.md):
  Plots scores of gene sets.

- [`plotQCvsAnnotation`](https://ccb-hms.github.io/scDiagnostics/reference/plotQCvsAnnotation.md):
  Plots QC metrics versus annotations.

## Misc

Miscellaneous functions for various tasks.

- [`processPCA`](https://ccb-hms.github.io/scDiagnostics/reference/processPCA.md):
  Process
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  objects to compute PCA.

- [`projectPCA`](https://ccb-hms.github.io/scDiagnostics/reference/projectPCA.md):
  Projects new data into PCA space.

- [`projectSIR`](https://ccb-hms.github.io/scDiagnostics/reference/projectSIR.md):
  Projects new data into SIR space.

- [`calculateCategorizationEntropy`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCategorizationEntropy.md):
  Calculates categorization entropy for clusters.

The package is built to facilitate in-depth analysis and visualization
of single-cell data, enhancing the understanding of cell type
similarities and differences across datasets.

## See also

Useful links:

- <https://github.com/ccb-hms/scDiagnostics>

- <https://ccb-hms.github.io/scDiagnostics>

- Report bugs at <https://github.com/ccb-hms/scDiagnostics/issues>

## Author

**Maintainer**: Anthony Christidis
<anthony-alexander_christidis@hms.harvard.edu>
([ORCID](https://orcid.org/0000-0002-4565-6279))

Authors:

- Andrew Ghazi

- Smriti Chawla

- Ludwig Geistlinger

- Robert Gentleman

Other contributors:

- Nitesh Turaga \[contributor\]
