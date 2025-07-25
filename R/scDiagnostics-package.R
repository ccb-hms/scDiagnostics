#' scDiagnostics: Single-Cell Diagnostics Package
#'
#' @description
#' `scDiagnostics` is a comprehensive toolkit designed for the analysis and diagnostics of single-cell RNA sequencing (scRNA-seq) data. This package provides functionalities for comparing principal components, visualizing canonical correlation analysis (CCA) outputs, and plotting cell type-specific MDS and PCA projections.
#'
#' @details
#' The package includes the following key functionalities, organized by their specific purposes:
#'
#' @section Visualization of Cell Type Annotations:
#' Functions for visualizing differences between query and reference datasets across multiple cell types.
#' \itemize{
#'   \item \code{\link{boxplotPCA}}: Boxplots of PCA scores for cell types.
#'   \item \code{\link{calculateDiscriminantSpace}}: Calculates discriminant space, with a plot method for visualization.
#'   \item \code{\link{plotCellTypeMDS}}: Creates MDS plots for cell types using query and reference datasets.
#'   \item \code{\link{plotCellTypePCA}}: Plots principal components for different cell types.
#' }
#'
#' @section Evaluation of Dataset Alignment:
#' Functions for visualizing differences between query and reference datasets for a specific cell type.
#' \itemize{
#'   \item \code{\link{calculateGraphIntegration}}: Calculate a graph network and form communities of cells.
#'   \item \code{\link{calculateWassersteinDistance}}: Wasserstein distance for different cell types.
#'   \item \code{\link{comparePCA}}: Compares PCA results, with a plot method for visualization.
#'   \item \code{\link{comparePCASubspace}}: Compares PCA subspace, with a plot method for visualization.
#'   \item \code{\link{plotPairwiseDistancesDensity}}: Plots the density of pairwise distances.
#' }
#'
#' @section Calculation of Statistical Measures to Compare Two Datasets:
#' Functions to compute statistical measures to compare two datasets.
#' \itemize{
#'   \item \code{\link{calculateAveragePairwiseCorrelation}}: Calculates average pairwise correlation, with a plot method for visualization.
#'   \item \code{\link{calculateCramerPValue}}: Calculates the p-value using Cramer's test.
#'   \item \code{\link{calculateHotellingPValue}}: Calculates the p-value using Hotelling's T-squared test.
#'   \item \code{\link{calculateMMDPValue}}: Calculates the p-value using maximum mean discrepancy.
#'   \item \code{\link{regressPC}}: Performs regression on principal components, with a plot method for visualization.
#' }
#'
#' @section Evaluation of Marker Gene Alignment:
#' Functions for calculating overlap measures of genes between two datasets.
#' \itemize{
#'   \item \code{\link{calculateHVGOverlap}}: Calculates overlap of highly variable genes (HVG) between datasets.
#'   \item \code{\link{calculateVarImpOverlap}}: Calculates overlap of variable importance measures between datasets.
#'   \item \code{\link{compareMarkers}}: Compares marker genes between datasets.
#' }
#'
#' @section Visualization of Marker Expressions:
#' Functions for to visualize and compare the expression of markers between a reference and a query dataset.
#' \itemize{
#'   \item \code{\link{plotGeneExpressionDimred}}: Plots gene expression in a dimensionality reduction space.
#'   \item \code{\link{plotMarkerExpression}}: Plots marker expression levels.
#' }
#'
#' @section Anomaly Detection (Global and Cell Type-Specific):
#' Functions for detecting anomalies at both the global and cell type-specific levels.
#' \itemize{
#'   \item \code{\link{detectAnomaly}}: Detects anomalies in the data, with a plot method for visualization.
#'   \item \code{\link{calculateCellSimilarityPCA}}: Calculates cell similarity in PCA space, with a plot method for visualization.
#' }
#'
#' @section Calculation of Distances Between Specific Cells and Cell Populations:
#' Functions for calculating distances between specific cells and cell populations.
#' \itemize{
#'   \item \code{\link{calculateCellDistances}}: Calculates distances between cells, with a plot method for visualization.
#'   \item \code{\link{calculateCellDistancesSimilarity}}: Calculates similarity based on cell distances, with a plot method for visualization.
#' }
#'
#' @section Visualization of QC and Annotation Scores:
#' Functions for visualizing quality control (QC) metrics or other characteristics of the data.
#' \itemize{
#'   \item \code{\link{histQCvsAnnotation}}: Plots histograms of QC metrics versus annotations.
#'   \item \code{\link{plotGeneSetScores}}: Plots scores of gene sets.
#'   \item \code{\link{plotQCvsAnnotation}}: Plots QC metrics versus annotations.
#' }
#'
#' @section Misc:
#' Miscellaneous functions for various tasks.
#' \itemize{
#'   \item \code{\link{processPCA}}: Process \code{\linkS4class{SingleCellExperiment}} objects to compute PCA.
#'   \item \code{\link{projectPCA}}: Projects new data into PCA space.
#'   \item \code{\link{projectSIR}}: Projects new data into SIR space.
#'   \item \code{\link{calculateCategorizationEntropy}}: Calculates categorization entropy for clusters.
#' }
#'
#' The package is built to facilitate in-depth analysis and visualization of single-cell data, enhancing the understanding of cell type similarities and differences across datasets.
#'
#' @name scDiagnostics-package
#' @aliases scDiagnostics
#' @keywords package internal
"_PACKAGE"
