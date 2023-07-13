#' Exploring Cell-Type Specific Patterns: Counts per Cell and Annotation Scores Visualization
#'
#' Creates histograms of gene expression values and scores obtained from cell type annotation methods for a specific cell type.
#'
#' @param cell_type_scores A numeric vector of scores for the specific cell type obtained from a cell type annotation method.
#' @param query_data An object of class "SingleCellExperiment" containing a numeric expression matrix.
#' @param cell_type A character string specifying the name of the cell type for which the distribution of log-transformed counts per cell is to be visualized.
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#'
#' @return A grid.arrange object displaying histograms of log-transformed counts and scores.
#'         This object can be further customized or used for additional plot manipulations.
#' @export
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
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # Log-transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR
#' cell_type_scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)$scores[, "CD4"]
#'
#' # Visualize the distribution of log-transformed counts and scores
#' plotCellTypeDistribution(cell_type_scores, query_data, "CD4")
#'
#' # Note: Users can use any cell type annotation method of their choice to obtain the scores.
#' # Ensure that the scores and log-transformed counts are provided to the function for visualization.
#'
plotCellTypeDistribution <- function(cell_type_scores, query_data, cell_type) {

  scores_df <- data.frame(Scores = cell_type_scores)

  # Pre-process logcounts
  log_expr <- as.matrix(assay(query_data, "logcounts"))
  counts_per_cell <- data.frame(CountsPerCell = colSums(log_expr))

  # Create a plot object for scores histogram
  scores_plot <- ggplot(scores_df, aes(x = Scores)) +
    geom_histogram(color = "black", fill = "white") +
    xlab("Scores") +
    ylab("Frequency") +
    theme_bw()

  # Create a plot object for total logcounts histogram
  counts_per_cell_plot <- ggplot(counts_per_cell, aes(x = CountsPerCell)) +
    geom_histogram(color = "black", fill = "white") +
    xlab("Total Logcounts") +
    ylab("Frequency") +
    theme_bw()

  return(grid.arrange(scores_plot, counts_per_cell_plot, ncol = 2))
}
