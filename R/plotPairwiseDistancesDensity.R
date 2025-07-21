#' @title Ridgeline Plot of Pairwise Distance Analysis
#'
#' @description
#' This function calculates pairwise distances or correlations between query and reference cells of a specified cell type
#' and visualizes the results using ridgeline plots, displaying the density distribution for each comparison.
#'
#' @details
#' Designed for \code{\linkS4class{SingleCellExperiment}} objects, this function subsets data for specified cell types,
#' computes pairwise distances or correlations, and visualizes these measurements through ridgeline plots.
#' The plots help evaluate the consistency and differentiation of annotated cell types within single-cell datasets.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} containing the single-cell
#' expression data and metadata.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing the single-cell
#' expression data and metadata.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_type_query The query cell type for which distances or correlations are calculated.
#' @param cell_type_ref The reference cell type for which distances or correlations are calculated.
#' @param pc_subset A numeric vector specifying which principal components to use in the analysis. Default is 1:5.
#' If set to \code{NULL}, the assay data is used directly for computations without dimensionality reduction.
#' @param distance_metric The distance metric to use for calculating pairwise distances, such as euclidean, manhattan, etc.
#'                        Set to "correlation" to calculate correlation coefficients.
#' @param correlation_method The correlation method to use when \code{distance_metric} is "correlation".
#'                           Possible values are "pearson" and "spearman".
#' @param bandwidth Numeric value controlling the smoothness of the density estimate; smaller values create more detailed curves. Default is 0.25.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A ggplot2 object showing ridgeline plots of calculated distances or correlations.
#'
#' @export
#'
#' @seealso \code{\link{calculateWassersteinDistance}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Example usage of the function
#' plotPairwiseDistancesDensity(query_data = query_data,
#'                              reference_data = reference_data,
#'                              query_cell_type_col = "SingleR_annotation",
#'                              ref_cell_type_col = "expert_annotation",
#'                              cell_type_query = "CD8",
#'                              cell_type_ref = "CD8",
#'                              pc_subset = 1:5,
#'                              distance_metric = "euclidean",
#'                              correlation_method = "pearson")
#'
#' @importFrom stats cor dist
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom ggridges geom_density_ridges
#'
# Plot density curves for distances between cells
plotPairwiseDistancesDensity <- function(
        query_data,
        reference_data,
        query_cell_type_col,
        ref_cell_type_col,
        cell_type_query,
        cell_type_ref,
        pc_subset = 1:5,
        distance_metric = c("correlation", "euclidean"),
        correlation_method = c("spearman", "pearson"),
        bandwidth = 0.25,
        assay_name = "logcounts",
        max_cells = 2500) {

    # Match argument for distance_metric
    distance_metric <- match.arg(distance_metric)

    # Match argument for correlation method
    correlation_method <- match.arg(correlation_method)

    # Check if bandwidth is valid
    if (!is.numeric(bandwidth) || length(bandwidth) != 1 || bandwidth <= 0 || bandwidth > 2) {
        stop("\'bandwidth\' must be a single positive numeric value between 0 and 2.")
    }

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = c(cell_type_query, cell_type_ref),
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)

    # Convert to matrix and potentially applied PCA dimensionality reduction
    if(!is.null(pc_subset)){
        # Project query data onto PCA space of reference data
        pca_output <- projectPCA(query_data = query_data,
                                 reference_data = reference_data,
                                 query_cell_type_col = query_cell_type_col,
                                 ref_cell_type_col = ref_cell_type_col,
                                 pc_subset = pc_subset,
                                 assay_name = assay_name)
        ref_mat <- pca_output[pca_output[["dataset"]] == "Reference" &
                                  pca_output[["cell_type"]] == cell_type_ref,
                              paste0("PC", pc_subset)]
        query_mat <- pca_output[pca_output[["dataset"]] == "Query" &
                                    pca_output[["cell_type"]] == cell_type_query,
                                paste0("PC", pc_subset)]
    } else{
        # Subset query data to the specified cell type
        query_data_subset <- query_data[, !is.na(query_data[[query_cell_type_col]]) &
                                            query_data[[query_cell_type_col]] == cell_type_query]
        query_mat <- t(as.matrix(assay(query_data_subset, assay_name)))
        ref_mat <- t(as.matrix(assay(reference_data, assay_name)))
    }

    # Combine query and reference matrices
    combined_mat <- rbind(query_mat, ref_mat)

    # Calculate pairwise distances or correlations for all comparisons
    if (distance_metric == "correlation") {
        dist_matrix <- cor(t(combined_mat), method = correlation_method)
    } else {
        dist_matrix <- as.matrix(dist(combined_mat, method = distance_metric))
    }

    # Extract the distances or correlations for the different pairwise comparisons
    num_query_cells <- nrow(query_mat)
    num_ref_cells <- nrow(ref_mat)
    query_indices <- seq_len(num_query_cells)
    ref_indices <- seq_len(num_ref_cells) + num_query_cells
    dist_query_query <- dist_matrix[query_indices, query_indices]
    dist_ref_ref <- dist_matrix[ref_indices, ref_indices]
    dist_query_ref <- dist_matrix[query_indices, ref_indices]

    # Create data frame for plotting
    ref_ref_name <- paste("Ref", cell_type_ref,
                          "vs", "Ref", cell_type_ref)
    query_query_name <- paste("Query", cell_type_query,
                              "vs", "Query", cell_type_query)
    query_ref_name <- paste("Query", cell_type_query,
                            "vs", "Ref", cell_type_ref)
    dist_df <- data.frame(
        Comparison = factor(c(
            rep(ref_ref_name, length(dist_ref_ref)),
            rep(query_query_name, length(dist_query_query)),
            rep(query_ref_name, length(dist_query_ref))),
            levels = c(ref_ref_name,
                       query_query_name,
                       query_ref_name)),
        Distance = c(as.vector(dist_ref_ref),
                     as.vector(dist_query_query),
                     as.vector(dist_query_ref))
    )

    # Plot ridge line plots
    ridgeline_plot <- ggplot2::ggplot(
        dist_df, ggplot2::aes(
            x = .data[["Distance"]], y = .data[["Comparison"]],
            fill = .data[["Comparison"]])) +
        ggridges::geom_density_ridges(scale = 1.5,
                                      bandwidth = bandwidth) +
        ggplot2::scale_fill_manual(values = setNames(
            c("#5DADE2", "#BA55D3", "#D9534F"),
            c(ref_ref_name, query_query_name, query_ref_name)
        )) +
        ggplot2::labs(
            x = ifelse(distance_metric == "correlation",
                       ifelse(correlation_method == "spearman",
                              "Spearman Correlation",
                              "Pearson Correlation"),
                       "Distance"),
            y = NULL,
            title = "Pairwise Distance Analysis") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold",
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            legend.position = "bottom") +
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE,
                                                     nrow = 3))

    # Adjust x-axis if correlation is used
    if(distance_metric == "correlation"){
        ridgeline_plot <- ridgeline_plot + ggplot2::xlim(-1, 1)
    }

    # Return the plot
    return(ridgeline_plot)
}
