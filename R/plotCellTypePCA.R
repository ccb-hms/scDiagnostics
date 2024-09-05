#' @title Plot Principal Components for Different Cell Types
#'
#' @description 
#' This function plots the principal components for different cell types in the query and reference datasets.
#'
#' @details
#' This function projects the query dataset onto the principal component space of the reference dataset and then plots the 
#' specified principal components for the specified cell types.
#' It uses the `projectPCA` function to perform the projection and \code{ggplot2} to create the plots.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is 1:5.
#'
#' @return A ggplot object representing the boxplots of specified principal components for the given cell types and datasets.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#' 
#' # Plot the PC data
#' pc_plot <- plotCellTypePCA(query_data = query_data, 
#'                            reference_data = reference_data,
#'                            cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                            query_cell_type_col = "expert_annotation", 
#'                            ref_cell_type_col = "expert_annotation", 
#'                            pc_subset = 1:5)
#' pc_plot
#' 
#' @importFrom stats approxfun cancor density setNames
#' @importFrom utils combn
#'                          
# Function to plot PC for different cell types
plotCellTypePCA <- function(query_data, 
                            reference_data, 
                            query_cell_type_col, 
                            ref_cell_type_col, 
                            cell_types = NULL,
                            pc_subset = 1:5){
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset)
    
    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }
    
    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data, 
                             reference_data = reference_data, 
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col,
                             pc_subset = pc_subset)
    pca_output <- na.omit(pca_output)

    # Create all possible pairs of specified PCs
    plot_names <- paste0("PC", pc_subset, " (", sprintf("%.1f%%", attributes(reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset]), ")")
    pairs <- expand.grid(x = plot_names, y = plot_names)
    pairs <- pairs[pairs[["x"]] != pairs[["y"]], ]
    # Create a new data frame with all possible pairs of specified PCs
    data_pairs_list <- lapply(seq_len(nrow(pairs)), function(i) {
        x_col <- pairs[["x"]][i]
        y_col <- pairs[["y"]][i]
        data_frame <- data.frame(pca_output[, c(pairs[["x"]][i], pairs[["y"]][i])], 
                                 paste(pca_output[["dataset"]], pca_output[["cell_type"]], sep = " "))
        colnames(data_frame) <- c("x_value", "y_value", "cell_type_dataset")
        data_frame[["x"]] <- x_col
        data_frame[["y"]] <- y_col
        return(data_frame)
    })
    # Plot data
    data_pairs <- do.call(rbind, data_pairs_list)
    # Remove redundant data (to avoid duplicated plots)
    data_pairs <- data_pairs[as.numeric(data_pairs[["x"]]) < as.numeric(data_pairs[["y"]]),]
    
    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)), rep(sort(cell_types), each = 2))
    data_pairs[["cell_type_dataset"]] <- factor(data_pairs[["cell_type_dataset"]], levels = order_combinations)
    cell_type_colors <- generateColors(order_combinations, paired = TRUE)

    # Create the ggplot object (with facets if PCA)
    plot_obj <- ggplot2::ggplot(
        data_pairs, ggplot2::aes(x = .data[["x_value"]], 
                                 y = .data[["y_value"]],
                                 color = .data[["cell_type_dataset"]])) +
        ggplot2::geom_point(alpha = 0.5, size = 1) +
        ggplot2::scale_color_manual(values = cell_type_colors, 
                                    name = "Cell Types") + 
        ggplot2::facet_grid(rows = ggplot2::vars(.data[["y"]]), 
                            cols = ggplot2::vars(.data[["x"]]), 
                            scales = "free") +
        ggplot2::xlab("") + ggplot2::ylab("") + 
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
    
    # Return the plot
    return(plot_obj)
}


