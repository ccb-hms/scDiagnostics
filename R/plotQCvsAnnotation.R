#' Scatter plot: QC stats vs Cell Type Annotation Scores
#'
#' Creates a scatter plot to visualize the relationship between QC stats (e.g., library size) 
#' and cell type annotation scores for one or more cell types.
#'
#' @details This function generates a scatter plot to explore the relationship between various quality 
#' control (QC) statistics, such as library size and mitochondrial percentage, and cell type 
#' annotation scores. By examining these relationships, users can assess whether specific QC 
#' metrics, systematically influence the confidence in cell type annotations, 
#' which is essential for ensuring reliable cell type annotation.
#' 
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} containing the single-cell 
#' expression data and metadata.
#' @param qc_col character. A column name in the \code{colData} of \code{query_data} that 
#' contains the QC stats of interest.
#' @param label_col character. The column name in the \code{colData} of \code{query_data} 
#' that contains the cell type labels.
#' @param score_col character. The column name in the \code{colData} of \code{query_data} that 
#' contains the cell type annotation scores.
#' @param label character. A vector of cell type labels to plot (e.g., c("T-cell", "B-cell")).  
#' Defaults to \code{NULL}, which will include all the cells.
#'
#' @return A ggplot object displaying a scatter plot of QC stats vs annotation scores, 
#'         where each point represents a cell, color-coded by its cell type.
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
#' set.seed(100)
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
#'plotQCvsAnnotation(query_data = query_data, 
#'                   qc_col = "percent.mito", 
#'                   label_col = "labels", 
#'                   score_col = "cell_scores", 
#'                   label = NULL)
#'                   
#' # Generate scatter plots
#' plotQCvsAnnotation(query_data = query_data, 
#'                    qc_col = "percent.mito", 
#'                    label_col = "labels", 
#'                    score_col = "cell_scores", 
#'                    label = c("CD4", "CD8"))
#'                    
#' @import ggplot2
#' @export
#'
plotQCvsAnnotation <- function(query_data, 
                               qc_col, 
                               label_col, 
                               score_col, 
                               label = NULL) {
  
  # Sanity checks
  
  # Check if query_data is a SingleCellExperiment object
  if (!is(query_data, "SingleCellExperiment")) {
    stop("query_data must be a SingleCellExperiment object.")
  }
  
  # Check if qc_col is a valid column name in query_data
  if (!qc_col %in% colnames(colData(query_data))) {
    stop("qc_col: '", qc_col, "' is not a valid column name in query_data.")
  }
  
  # Check if label_col is a valid column name in query_data
  if (!label_col %in% colnames(colData(query_data))) {
    stop("label_col: '", label_col, "' is not a valid column name in query_data.")
  }
  
  # Check if score_col is a valid column name in query_data
  if (!score_col %in% colnames(colData(query_data))) {
    stop("score_col: '", score_col, "' is not a valid column name in query_data.")
  }
  
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
  data <- data.frame(QCStats = qc_stats, 
                     Scores = cell_type_scores, 
                     CellType = cell_labels)
  
  # Create a scatter plot with color-coded points based on cell types or labels
  plot <- ggplot(data, aes(x = QCStats, 
                           y = Scores, 
                           color = CellType)) +
    geom_point() +
    xlab("QC stats") +
    ylab("Annotation Scores") +
    theme_bw()
  
  return(plot)
}