#' @title Histograms: QC Stats and Annotation Scores Visualization
#'
#' @description
#' This function generates histograms for visualizing the distribution of quality control (QC) statistics and 
#' annotation scores associated with cell types in single-cell genomic data. 
#' 
#' @details The particularly useful in the analysis of data from single-cell experiments, 
#' where understanding the distribution of these metrics is crucial for quality assessment and 
#' interpretation of cell type annotations.
#'
#' @param query_data  A \code{\linkS4class{SingleCellExperiment}} containing the single-cell 
#' expression data and metadata.
#' @param qc_col character. A column name in the \code{colData} of \code{query_data} that 
#' contains the QC stats of interest.
#' @param label_col character. The column name in the \code{colData} of \code{query_data} 
#' that contains the cell type labels.
#' @param score_col numeric. The column name in the \code{colData} of \code{query_data} that 
#' contains the cell type scores.
#' @param label character. A vector of cell type labels to plot (e.g., c("T-cell", "B-cell")).  
#' Defaults to \code{NULL}, which will include all the cells.
#'
#' @return A object containing two histograms displayed side by side. 
#' The first histogram represents the distribution of QC stats, 
#' and the second histogram represents the distribution of annotation scores.
#' 
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(SingleR)
#' library(gridExtra)
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
#' 
#' # Assign labels to query data
#' colData(query_data)$labels <- pred$labels
#' 
#' # Get annotation scores
#' scores <- apply(pred$scores, 1, max)
#'
#' # Assign scores to query data
#' colData(query_data)$cell_scores <- scores
#'
#' # Generate histograms
#' histQCvsAnnotation(query_data = query_data, 
#'                   qc_col = "percent.mito", 
#'                   label_col = "labels", 
#'                   score_col = "cell_scores", 
#'                   label = c("CD4", "CD8"))
#'                   
#' histQCvsAnnotation(query_data = query_data, 
#'                    qc_col = "percent.mito", 
#'                    label_col = "labels", 
#'                    score_col = "cell_scores", 
#'                    label = NULL)
#'
#' @export
histQCvsAnnotation <- function(query_data, 
                               qc_col = qc_col, 
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

  # Combine QC stats, scores, and labels into a data frame
  data <- data.frame(QCStats = qc_stats, Scores = cell_type_scores)
  
  # Create histogram for QC stats
  qc_histogram <- ggplot2::ggplot(data, aes(x = QCStats)) +
      ggplot2::geom_histogram(color = "black", fill = "white") +
      ggplot2::xlab(paste(qc_col)) +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw()
  
  # Create histogram for scores
  scores_histogram <- ggplot2::ggplot(data, aes(x = Scores)) +
      ggplot2::geom_histogram(color = "black", fill = "white") +
      ggplot2::xlab("Annotation Scores") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw()
  
  # Return the list of plots
  return(gridExtra::grid.arrange(qc_histogram, scores_histogram, ncol = 2))
}
