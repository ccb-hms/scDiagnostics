#' Scatter plot: Total UMIs vs. Annotation Scores
#'
#' Creates a scatter plot to visualize the relationship between total UMIs (Unique Molecular Identifiers) and annotation scores for a specific cell type.
#'
#' @param query_data An object of class "SingleCellExperiment" containing the single-cell expression data and metadata.
#' @param cell_type_scores A numeric vector of scores for the specific cell type obtained from a cell type annotation method.
#' @param cell_type_labels A character string representing the column name in the colData(query_data) that contains the cell type labels.
#' @param cell_type A character string specifying the name of the cell type for which the scatter plot is to be generated.
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
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Extract scores and labels for a specific cell type
#' cell_type <- "CD4"
#' score <- pred$scores[pred$labels == cell_type, cell_type]
#'
#' # Assign labels to query data
#' colData(query_data)$labels <- scores$labels
#'
#' # Generate scatter plot for Total UMIs vs. Annotation Scores
#' plotUMIsAnnotationScatter(query_data, score, "labels", "CD4")
#'
#' # Note: Users can use any cell type annotation method of their choice to obtain the cell type scores.
#' 
#'
plotUMIsAnnotationScatter <- function(query_data, cell_type_scores, cell_type_labels, cell_type) {

  # Pre-process expression matrix and extract scores for a specific cell type
  indx <- query_data[[cell_type_labels]] == cell_type
  total_umis <- colSums(as.matrix(assay(query_data, "logcounts")[,indx]))

  # Combine total UMIs and annotation scores into a data frame
  data <- data.frame(TotalUMIs = total_umis, Scores = cell_type_scores)

  # Create a scatter plot
  plot <- ggplot(data, aes(x = TotalUMIs, y = Scores)) +
    geom_point() +
    xlab("Total logcounts") +
    ylab("Annotation Scores") +
    theme_bw()

  return(plot)
}
