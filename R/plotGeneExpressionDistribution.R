#' Plot gene expression distribution from overall and cell type-specific perspective
#'
#' This function generates histogram plots to visualize the distribution of gene expression values for a specific gene
#' across the overall dataset and within a specified cell type.
#'
#' @param se_object An object of class "SingleCellExperiment" containing log-transformed expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param cell_type_labels A character column name from the colData table of the SingleCellExperiment object representing cell type information.
#' @param cell_type A character string representing the name of the cell type for which the specific gene distribution is to be visualized.
#' @param feature A character string representing the gene name for which the distribution is to be visualized.
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @import SingleCellExperiment
#'
#' @return A ggplot object displaying the histogram plots of gene expression distribution from the overall dataset and cell type-specific perspective.
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
#' # Log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR or any other method
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#'
#' # # Visualize the distribution of log-transformed counts and scores for the se_object
#' # Note: Users can use SingleR or any other method to obtain the cell type scores.
#' # Ensure that the scores and log-transformed counts are provided to the function for visualization.
#' # For demonstration we have used query_data.
#' plotGeneExpressionDistribution(query_data, "labels", "B_and_plasma", "VPREB3")
#'
plotGeneExpressionDistribution <- function(query_data, cell_type_labels, cell_type, feature) {
  # Get expression of the specified feature
  overall <- as.data.frame(assay(query_data, "logcounts")[feature, ])
  indx <- which(colData(query_data)[, cell_type_labels] == cell_type)
  cell_specific <- as.data.frame(assay(query_data, "logcounts")[feature, indx])

  # Assign a common column name to expression values
  colnames(overall) <- colnames(cell_specific) <- feature

  # Pre-process the data
  overall$group <- "Overall Distribution"
  cell_specific$group <- "Cell type-specific Distribution"

  # Combine overall and cell-specific data frames
  df <- rbind(overall, cell_specific)

  # Generate histograms for overall distribution and cell type-specific distribution of the gene
  plots <- ggplot(df, aes_string(x = feature)) +
    geom_histogram(binwidth = 0.5, color = "black", fill = "white") +
    facet_wrap(
      ~ group,
      scales = "free_y"
    ) +
    labs(
      title = cell_type,
      x = paste("Log gene expression", feature),
      y = "Frequency"
    ) +
    theme_bw()

  return(plots)
}
