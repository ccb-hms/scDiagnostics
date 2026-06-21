# Package index

## Visualization of Cell Type Annotations

- [`boxplotPCA()`](https://ccb-hms.github.io/scDiagnostics/reference/boxplotPCA.md)
  : Plot Principal Components for Different Cell Types
- [`calculateDiscriminantSpace()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateDiscriminantSpace.md)
  [`plot(`*`<calculateDiscriminantSpaceObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateDiscriminantSpace.md)
  : Project Query Data onto a Unified Discriminant Space of Reference
  Data
- [`calculateSIRSpace()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateSIRSpace.md)
  [`plot(`*`<calculateSIRSpaceObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateSIRSpace.md)
  : Calculate Sliced Inverse Regression (SIR) Space for Different Cell
  Types
- [`plotCellTypeMDS()`](https://ccb-hms.github.io/scDiagnostics/reference/plotCellTypeMDS.md)
  : Plot Reference and Query Cell Types using MDS
- [`plotCellTypePCA()`](https://ccb-hms.github.io/scDiagnostics/reference/plotCellTypePCA.md)
  : Plot Principal Components for Different Cell Types

## Evaluation of Dataset Alignment

- [`calculateGraphIntegration()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateGraphIntegration.md)
  [`plot(`*`<calculateGraphIntegrationObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateGraphIntegration.md)
  : Calculate Graph Community Integration Diagnostics
- [`calculateWassersteinDistance()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateWassersteinDistance.md)
  [`plot(`*`<calculateWassersteinDistanceObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateWassersteinDistance.md)
  : Compute Wasserstein Distance Distributions Between Query and
  Reference Datasets
- [`comparePCA()`](https://ccb-hms.github.io/scDiagnostics/reference/comparePCA.md)
  [`plot(`*`<comparePCAObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/comparePCA.md)
  : Compare Principal Components Analysis (PCA) Results
- [`comparePCASubspace()`](https://ccb-hms.github.io/scDiagnostics/reference/comparePCASubspace.md)
  [`plot(`*`<comparePCASubspaceObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/comparePCASubspace.md)
  : Compare Subspaces Spanned by Top Principal Components
- [`plotPairwiseDistancesDensity()`](https://ccb-hms.github.io/scDiagnostics/reference/plotPairwiseDistancesDensity.md)
  : Ridgeline Plot of Pairwise Distance Analysis

## Calculation of Statistical Measures to Compare Two Datasets

- [`calculateAveragePairwiseCorrelation()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateAveragePairwiseCorrelation.md)
  [`plot(`*`<calculateAveragePairwiseCorrelationObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateAveragePairwiseCorrelation.md)
  : Compute Average Pairwise Correlation between Cell Types
- [`calculateCramerPValue()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCramerPValue.md)
  : Calculate Cramer Test P-Values for Two-Sample Comparison of
  Multivariate ECDFs
- [`calculateHotellingPValue()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateHotellingPValue.md)
  : Perform Hotelling's T-squared Test on PCA Scores for Single-cell
  RNA-seq Data
- [`calculateMMDPValue()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateMMDPValue.md)
  : Calculate Maximum Mean Discrepancy P-Values for Two-Sample
  Comparison
- [`plot(`*`<regressPCObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/regressPC.md)
  [`regressPC()`](https://ccb-hms.github.io/scDiagnostics/reference/regressPC.md)
  : Plot Regression Results on Principal Components

## Anomaly Detection (Global and Cell Type-Specific)

- [`detectAnomaly()`](https://ccb-hms.github.io/scDiagnostics/reference/detectAnomaly.md)
  [`plot(`*`<detectAnomalyObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/detectAnomaly.md)
  : PCA Anomaly Scores via Isolation Forests with Visualization
- [`calculateReconstructionError()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateReconstructionError.md)
  [`plot(`*`<calculateReconstructionErrorObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateReconstructionError.md)
  : Calculate PCA Reconstruction Errors for Out-of-Distribution Anomaly
  Detection
- [`calculateCellSimilarityPCA()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellSimilarityPCA.md)
  [`plot(`*`<calculateCellSimilarityPCAObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellSimilarityPCA.md)
  : Calculate Cell Similarity Using PCA Loadings

## Calculation of Distances Between Specific Cells and Cell Populations

- [`calculateCellDistances()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellDistances.md)
  [`plot(`*`<calculateCellDistancesObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellDistances.md)
  : Compute Cell Distances Between Reference and Query Data
- [`calculateCellDistancesSimilarity()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCellDistancesSimilarity.md)
  : Function to Calculate Bhattacharyya Coefficients and Hellinger
  Distances

## Evaluation of Marker Gene Alignment

- [`calculateHVGOverlap()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateHVGOverlap.md)
  : Calculate the Overlap Coefficient for Highly Variable Genes
- [`calculateGeneShifts()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateGeneShifts.md)
  [`plot(`*`<calculateGeneShiftsObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/calculateGeneShifts.md)
  : Calculate Top Loading Gene Expression Shifts
- [`calculateVarImpOverlap()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateVarImpOverlap.md)
  : Compare Gene Importance Across Datasets Using Random Forest
- [`compareMarkers()`](https://ccb-hms.github.io/scDiagnostics/reference/compareMarkers.md)
  [`plot(`*`<compareMarkersObject>`*`)`](https://ccb-hms.github.io/scDiagnostics/reference/compareMarkers.md)
  : Compare Marker Gene Expression between Query and Reference Data

## Visualization of Marker Expressions

- [`plotGeneExpressionDimred()`](https://ccb-hms.github.io/scDiagnostics/reference/plotGeneExpressionDimred.md)
  : Visualize gene expression on a dimensional reduction plot
- [`plotMarkerExpression()`](https://ccb-hms.github.io/scDiagnostics/reference/plotMarkerExpression.md)
  : Plot gene expression distribution from overall and cell
  type-specific perspective

## Visualization of QC and Annotation Scores

- [`histQCvsAnnotation()`](https://ccb-hms.github.io/scDiagnostics/reference/histQCvsAnnotation.md)
  : Histograms: QC Stats and Annotation Scores Visualization
- [`plotQCvsAnnotation()`](https://ccb-hms.github.io/scDiagnostics/reference/plotQCvsAnnotation.md)
  : Scatter plot: QC stats vs Cell Type Annotation Scores
- [`plotGeneSetScores()`](https://ccb-hms.github.io/scDiagnostics/reference/plotGeneSetScores.md)
  : Visualization of gene sets or pathway scores on dimensional
  reduction plot

## Misc

- [`processPCA()`](https://ccb-hms.github.io/scDiagnostics/reference/processPCA.md)
  : Process PCA for SingleCellExperiment Objects
- [`projectPCA()`](https://ccb-hms.github.io/scDiagnostics/reference/projectPCA.md)
  : Project Query Data Onto PCA Space of Reference Data
- [`projectSIR()`](https://ccb-hms.github.io/scDiagnostics/reference/projectSIR.md)
  : Project Query Data Onto SIR Space of Reference Data
- [`calculateCategorizationEntropy()`](https://ccb-hms.github.io/scDiagnostics/reference/calculateCategorizationEntropy.md)
  : Calculate Categorization Entropy

## Datasets

- [`reference_data`](https://ccb-hms.github.io/scDiagnostics/reference/reference_data.md)
  : Reference Single-Cell RNA-Seq Dataset
- [`query_data`](https://ccb-hms.github.io/scDiagnostics/reference/query_data.md)
  : Query Single-Cell RNA-Seq Dataset
- [`qc_data`](https://ccb-hms.github.io/scDiagnostics/reference/qc_data.md)
  : Quality Control Single-Cell RNA-Seq Dataset
