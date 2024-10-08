#' @title Pairwise Distance Analysis and Density Visualization
#'
#' @description
#' Calculates pairwise distances or correlations between query and reference cells of a specific cell type.
#' 
#' @details  
#' The function works with \code{\linkS4class{SingleCellExperiment}} objects, ensuring 
#' compatibility with common single-cell analysis workflows. It subsets the data for specified cell types, 
#' computes pairwise distances or correlations, and visualizes these measurements using density plots. By comparing the distances and correlations, 
#' one can evaluate the consistency and reliability of annotated cell types within single-cell datasets.
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
#' If set to \code{NULL} then no dimensionality reduction is performed and the assay data is used directly for computations.
#' @param distance_metric The distance metric to use for calculating pairwise distances, such as euclidean, manhattan etc.
#'                        Set it to "correlation" for calculating correlation coefficients.
#' @param correlation_method The correlation method to use when distance_metric is "correlation".
#'                           Possible values: "pearson", "spearman".
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#'
#' @return A plot generated by \code{ggplot2}, showing the density distribution of 
#'         calculated distances or correlations.
#'         
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
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
#' @export
#' 
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
        assay_name = "logcounts") {
    
    # Match argument for distance_metric
    distance_metric <- match.arg(distance_metric)
    
    # Match argument for correlation method
    correlation_method <- match.arg(correlation_method)
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = c(cell_type_query, cell_type_ref),
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)
    
    # Convert to matrix and potentially applied PCA dimensionality reduction
    if(!is.null(pc_subset)){
        # Project query data onto PCA space of reference data
        pca_output <- projectPCA(query_data = query_data, 
                                 reference_data = reference_data, 
                                 query_cell_type_col = query_cell_type_col,
                                 ref_cell_type_col = ref_cell_type_col,
                                 pc_subset = pc_subset,
                                 assay_name = assay_name)
        ref_mat <- pca_output[which(pca_output[["dataset"]] == "Reference" &
                                        pca_output[["cell_type"]] == cell_type_ref), 
                              paste0("PC", pc_subset)]
        query_mat <- pca_output[which(pca_output[["dataset"]] == "Query" &
                                          pca_output[["cell_type"]] == cell_type_query), 
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
        if (correlation_method == "pearson") {
            dist_matrix <- cor(t(combined_mat), method = "pearson")
        } else if (correlation_method == "spearman") {
            dist_matrix <- cor(t(combined_mat), method = "spearman")
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
    query_indices <- seq_len(num_query_cells)
    ref_indices <- seq_len(num_ref_cells) + num_query_cells
    dist_query_query <- dist_matrix[query_indices, query_indices]
    dist_ref_ref <- dist_matrix[ref_indices, ref_indices]
    dist_query_ref <- dist_matrix[query_indices, ref_indices]
    
    # Create data frame for plotting
    dist_df <- data.frame(
        Comparison = c(rep("Query vs Query", length(dist_query_query)),
                       rep("Reference vs Reference", length(dist_ref_ref)),
                       rep("Query vs Reference", length(dist_query_ref))),
        Distance = c(as.vector(dist_query_query),
                     as.vector(dist_ref_ref),
                     as.vector(dist_query_ref))
    )
    
    # Plot density plots with improved aesthetics
    ggplot2::ggplot(dist_df, ggplot2::aes(x = .data[["Distance"]], 
                                          color = .data[["Comparison"]])) +
        ggplot2::geom_density(alpha = 0.5, linewidth = 0.5, adjust = 2) +  
        ggplot2::scale_color_manual(values = c("#D9534F", "#BA55D3", "#5DADE2")) +
        ggplot2::labs(
            x = ifelse(distance_metric == "correlation", 
                       ifelse(correlation_method == "spearman", 
                              "Spearman Correlation", 
                              "Pearson Correlation"), 
                       "Distance"), 
            y = "Density", 
            title = "Pairwise Distance Analysis and Density Visualization") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray", 
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14, 
                                               face = "bold", 
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12), 
            axis.text = ggplot2::element_text(size = 10))
}
