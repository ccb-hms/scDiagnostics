#' Pairwise Distance Analysis and Density Visualization
#'
#' This function calculates pairwise distances or correlations between query and reference cells of a specific cell type based on the provided distance or correlation metric.
#' It then plots density plots to visualize the distribution of distances or correlations for different pairwise comparisons.
#'
#' @param query_data A SingleCellExperiment object containing numeric expression matrix for the query cells.
#' @param ref_data A SingleCellExperiment object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the query_data metadata specifying the cell types.
#' @param ref_cell_type_col The column name in the ref_data metadata specifying the cell types.
#' @param cell_type The cell type for which distances or correlations are calculated.
#' @param distance_metric The distance metric to use for calculating pairwise distances, such as "euclidean" or "manhattan" etc.
#'                        Set it to "correlation" for calculating correlation coefficients.
#' @param correlation_method The correlation method to use when distance_metric is "correlation".
#'                           Possible values: "pearson", "spearman".
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @import SingleCellExperiment
#' @export
#'
#' @examples
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
#' # log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#'
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- getTopHVGs(ref_data, n = 2000)
#' query_var <- getTopHVGs(query_data, n = 2000)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#'
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#'
#' # Example usage of the function
#' calculatePairwiseDistancesAndPlotDensity(query_data_subset, ref_data_subset, "labels", "reclustered.broad", "CD4", "euclidean")
#'
calculatePairwiseDistancesAndPlotDensity <- function(query_data, ref_data, query_cell_type_col, ref_cell_type_col, cell_type, distance_metric, correlation_method = "pearson") {
  # Subset query and reference data to the specified cell type
  query_data_subset <- query_data[, !is.na(query_data[[query_cell_type_col]]) & query_data[[query_cell_type_col]] == cell_type]
  ref_data_subset <- ref_data[, !is.na(ref_data[[ref_cell_type_col]]) & ref_data[[ref_cell_type_col]] == cell_type]

  # Convert to matrix
  query_mat <- t(as.matrix(assay(query_data_subset)))
  ref_mat <- t(as.matrix(assay(ref_data_subset)))

  # Combine query and reference matrices
  combined_mat <- rbind(query_mat, ref_mat)

  # Calculate pairwise distances or correlations for all comparisons
  if (distance_metric == "correlation") {
    if (correlation_method == "pearson") {
      dist_matrix <- cor(combined_mat, method = "pearson")
    } else if (correlation_method == "spearman") {
      dist_matrix <- cor(combined_mat, method = "spearman")
    } else {
      stop("Invalid correlation method. Available options: 'pearson', 'spearman'")
    }
  } else {
    dist_matrix <- dist(combined_mat, method = distance_metric)
  }

  # Convert dist_matrix to a square matrix
  dist_matrix <- as.matrix(dist_matrix)

  # Extract the distances or correlations for the different pairwise comparisons
  num_query_cells <- nrow(query_mat)
  num_ref_cells <- nrow(ref_mat)
  dist_query_query <- dist_matrix[1:num_query_cells, 1:num_query_cells]
  dist_ref_ref <- dist_matrix[(num_query_cells+1):(num_query_cells+num_ref_cells), (num_query_cells+1):(num_query_cells+num_ref_cells)]
  dist_query_ref <- dist_matrix[1:num_query_cells, (num_query_cells+1):(num_query_cells+num_ref_cells)]

  # Create data frame for plotting
  dist_df <- data.frame(
    Comparison = c(rep("Query vs Query", length(dist_query_query)),
                   rep("Reference vs Reference", length(dist_ref_ref)),
                   rep("Query vs Reference", length(dist_query_ref))),
    Distance = c(as.vector(dist_query_query),
                 as.vector(dist_ref_ref),
                 as.vector(dist_query_ref))
  )

  # Plot density plots
  ggplot(dist_df, aes(x = Distance, color = Comparison)) +
    geom_density() +
    labs(x = ifelse(distance_metric == "correlation", paste(correlation_method, "correlation"), "Distance"), y = "Density", title = "Pairwise Distance Analysis and Density Visualization") +
    theme_bw()
}
