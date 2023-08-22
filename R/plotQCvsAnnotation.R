#' Scatter plot: QC stats and Cell Type Scores
#'
#' Creates a scatter plot to visualize the relationship between QC stats (e.g., library size) and cell type scores for one or more cell types.
#'
#' @param query_data An object of class "SingleCellExperiment" containing the single-cell expression data and metadata.
#' @param qc_col A character string specifying the column name in the colData(query_data) that contains the QC stats of interest.
#' @param label_col A character string specifying the column name in the colData(query_data) that contains the cell type labels.
#' @param score_col A character string specifying the column name in the colData(query_data) that contains the cell type scores.
#' @param label A character vector of cell type labels to plot (e.g., c("T-cell", "B-cell")). If NULL, all cells will be included.
#'
#' @return A ggplot object displaying the scatter plot of total UMIs and annotation scores.
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(SingleR)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # Log-transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR
#' pred <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' pred <- as.data.frame(pred)
#'
#' # Assign labels to query data
#' colData(query_data)$labels <- pred$labels
#' scores <- apply(pred[,1:4], 1, max)
#'
#' # Assign scores to query data
#' colData(query_data)$cell_scores <- scores
#'
#' # Generate scatter plots
#' plotQCvsAnnotation(query_data = query_data, qc_col = "percent.mito", label_col = "labels", score_col = "cell_scores", label = c("CD4", "CD8"))
#' plotQCvsAnnotation(query_data = query_data, qc_col = "percent.mito", label_col = "labels", score_col = "cell_scores", label = NULL)
#'
plotQCvsAnnotation <- function(query_data, qc_col, label_col, score_col, label = NULL) {
  
  # Filter cells based on label if specified
  if (!is.null(label)) {
    index <- which(colData(query_data)[[label_col]] %in% label)
    query_data <- query_data[, index]
  }
  
  # Extract QC stats, scores, and labels
  qc_stats <- colData(query_data)[, qc_col]
  cell_type_scores <- colData(query_data)[, score_col]
  cell_labels <- colData(query_data)[[label_col]]
  
  # Combine QC stats, scores, and labels into a data frame
  data <- data.frame(QCStats = qc_stats, Scores = cell_type_scores, CellType = cell_labels)
  
  # Create a scatter plot with color-coded points based on cell types or labels
  plot <- ggplot(data, aes(x = QCStats, y = Scores, color = CellType)) +
    geom_point() +
    xlab("QC stats") +
    ylab("Annotation Scores") +
    theme_bw()
  
  return(plot)
}