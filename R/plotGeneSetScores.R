#' Visualization of gene sets or pathway scores on dimensional reduction plot
#'
#' Plot gene sets or pathway scores on PCA, TSNE, or UMAP. Single cells are color-coded by scores of gene sets or pathways.
#'
#' @param se_object An object of class "SingleCellExperiment" containing numeric expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param method A character string indicating the method for visualization ("PCA", "TSNE", or "UMAP").
#' @param feature A character string representing the name of the feature (score) in the colData(query data or reference data) to plot.
#'
#' @import scater
#' @importFrom scater plotPCA plotTSNE plotUMAP
#'
#' @return A plot object with the scores plotted on the selected dimensional reduction plot.
#' @export
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(AUCell)
#'
#' ## load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' ## log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Run PCA on the query data
#' query_data <- runPCA(query_data)
#'
#' # Assign scores to colData (users should ensure that the scores are present in the colData)
#' colData(query_data)$geneSetScores <- # ... (assign the scores here)
#'
#' # Plot gene set scores on PCA
#' plotGeneSetScores(query_data, method = "PCA", feature = "geneSetScores")
#'
#' # Note: Users can provide their own gene set scores in the colData of the 'se_object' object, using any method of their choice.
#' # Ensure that the scores are assigned to the colData and specify the correct feature name for visualization.
#' The example code shows the assignment of gene set scores to the colData of the query object for demonstration purposes.
#'
plotGeneSetScores <- function(se_object, method, feature) {

  # Check if the specified method is valid
  valid_methods <- c("PCA", "TSNE", "UMAP")
  if (!(method %in% valid_methods)) {
    stop("Invalid method. Please choose one of: ", paste(valid_methods, collapse = ", "))
  }

  # Create the plot object
  if (method == "PCA") {
    plot_obj <- plotPCA(query_data, colour_by = feature)
  } else if (method == "TSNE") {
    plot_obj <- plotTSNE(query_data, colour_by = feature)
  } else if (method == "UMAP") {
    plot_obj <- plotUMAP(query_data, colour_by = feature)
  }

  return(plot_obj)
}
