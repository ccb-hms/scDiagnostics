#' @title Visualize Principal Components for Different Cell Types
#'
#' @description 
#' This function plots the principal components for different cell types in the query and reference datasets.
#'
#' @details
#' This function projects the query dataset onto the principal component space of the reference dataset and then visualizes the 
#' specified principal components for the specified cell types.
#' It uses the `projectPCA` function to perform the projection and `ggplot2` to create the plots.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param n_components An integer specifying the number of principal components to use for projection. Defaults to 10. 
#' Must be less than or equal to the number of components available in the reference PCA.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#'
#' @return A ggplot object representing the boxplots of specified principal components for the given cell types and datasets.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load required libraries
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scran)
#' library(scater)
#'
#' # Load data (replace with your data loading)
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#' 
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#' 
#' # log transform datasets
#' ref_data <- scuttle::logNormCounts(ref_data)
#' query_data <- scuttle::logNormCounts(query_data)
#' 
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR::SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' 
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- scran::getTopHVGs(ref_data, n = 2000)
#' query_var <- scran::getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#'
#' # Run PCA on the reference data (assumed to be prepared)
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' pc_plot <- visualizeCellTypePCA(query_data_subset, ref_data_subset,
#'                                 n_components = 10,
#'                                 cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                                 query_cell_type_col = "labels", 
#'                                 ref_cell_type_col = "reclustered.broad", 
#'                                 pc_subset = c(1:5))
#' pc_plot
#' 
#' 
#' @importFrom stats approxfun cancor density setNames
#' @importFrom utils combn
#'                          
# Function to plot PC for different cell types
visualizeCellTypePCA <- function(query_data, reference_data, 
                                 n_components = 10, 
                                 cell_types = NULL,
                                 query_cell_type_col, 
                                 ref_cell_type_col, 
                                 pc_subset = c(1:5)){
    
    # Cell types
    if(is.null(cell_types)){
        cell_types <- na.omit(intersect(unique(query_data[[query_cell_type_col]]), unique(reference_data[[ref_cell_type_col]])))
    }
    
    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data, reference_data = reference_data, 
                             n_components = n_components, 
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col)
    pca_output <- na.omit(pca_output)

    # Create all possible pairs of specified PCs
    plot_names <- paste0("PC", pc_subset)
    pairs <- expand.grid(x = plot_names, y = plot_names)
    pairs <- pairs[pairs$x != pairs$y, ]
    # Create a new data frame with all possible pairs of specified PCs
    data_pairs_list <- lapply(1:nrow(pairs), function(i) {
        x_col <- pairs$x[i]
        y_col <- pairs$y[i]
        data_frame <- data.frame(pca_output[, c(x_col, y_col)], paste(pca_output$dataset, pca_output$cell_type, sep = " "))
        colnames(data_frame) <- c("x_value", "y_value", "cell_type_dataset")
        data_frame$x <- x_col
        data_frame$y <- y_col
        data_frame
    })
    # Plot data
    data_pairs <- do.call(rbind, data_pairs_list)
    # Remove redundant data (to avoid duplicated plots)
    data_pairs <- data_pairs[as.numeric(data_pairs$x) < as.numeric(data_pairs$y),]
    
    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)), rep(sort(cell_types), each = 2))
    data_pairs$cell_type_dataset <- factor(data_pairs$cell_type_dataset, levels = order_combinations)
    color_mapping <- setNames(RColorBrewer::brewer.pal(length(order_combinations), "Paired"), order_combinations)
    cell_type_colors <- color_mapping[order_combinations]

    # Create the ggplot object (with facets if PCA)
    plot_obj <- ggplot2::ggplot(data_pairs, ggplot2::aes(x = x_value, y = y_value, color = cell_type_dataset)) +
        ggplot2::geom_point(alpha = 0.5, size = 1) +
        ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Types") + 
        ggplot2::facet_grid(rows = ggplot2::vars(y), cols = ggplot2::vars(x), scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10))
    
    # Return the plot
    return(plot_obj)
}


