#' Histograms: QC Stats and Annotation Scores Visualization
#'
#' Creates histograms to visualize the distribution of QC stats and scores obtained from cell type annotation methods for one or more cell types.
#'
#' @param query_data An object of class "SingleCellExperiment" containing the single-cell expression data and metadata.
#' @param qc_col A character string specifying the column name in the colData(query_data) that contains the QC stats of interest.
#' @param label_col A character string specifying the column name in the colData(query_data) that contains the cell type labels.
#' @param score_col A character string specifying the column name in the colData(query_data) that contains the cell type scores.
#' @param label A character vector of cell type labels to plot (e.g., c("T-cell", "B-cell")). If NULL, all cells will be included.
#'
#' @import ggplot2
#'
#' #' @return A list of ggplot objects containing histograms of QC stats and annotation scores.
#' The list contains two ggplot objects: the first one displays the histogram of QC stats,
#' and the second one displays the histogram of annotation scores.
#'
#' @export
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
#' # Generate histograms
#' histQCvsAnnotation(query_data, "percent.mito", "labels", "cell_scores", c("CD4", "CD8"))
#' histQCvsAnnotation(query_data, "percent.mito", "labels", "cell_scores", NULL)
#'
histQCvsAnnotation <- function(query_data, qc_col, label_col, score_col, label = NULL) {
  
  # Filter cells based on label if specified
  if (!is.null(label)) {
    index <- which(colData(query_data)[[label_col]] %in% label)
    query_data <- query_data[, index]
  }
  
  # Extract QC stats, scores, and labels
  qc_stats <- colData(query_data)[, qc_col]
  cell_type_scores <- colData(query_data)[, score_col]
  cell_type_labels <- colData(query_data)[, label_col]
  
  # Combine QC stats, scores, and labels into a data frame
  data <- data.frame(QCStats = qc_stats, Scores = cell_type_scores)
  
  # Create histogram for QC stats
  qc_histogram <- ggplot(data, aes(x = QCStats)) +
    geom_histogram(color = "black", fill = "white") +
    xlab(paste(qc_col)) +
    ylab("Frequency") +
    theme_bw()
  
  # Create histogram for scores
  scores_histogram <- ggplot(data, aes(x = Scores)) +
    geom_histogram(color = "black", fill = "white") +
    xlab("Annotation Scores") +
    ylab("Frequency") +
    theme_bw()
  
  # Return a grid.arrange object displaying histograms of QC stats and annotation scores
  return(list(qc_histogram, scores_histogram))
}
