#' @title Create Faceted Scatter Plots for Specified PC Combinations From \code{detectAnomaly} Object
#'
#' @description This function generates faceted scatter plots for specified principal component (PC) combinations
#' within an anomaly detection object. It allows visualization of the relationship between specified
#' PCs and highlights anomalies detected by the Isolation Forest algorithm.
#'
#' @details The function extracts the specified PCs from the given anomaly detection object and generates
#' scatter plots for each pair of PCs. It uses \code{ggplot2} to create a faceted plot where each facet represents
#' a pair of PCs. Anomalies are highlighted in red, while normal points are shown in black.
#'
#' @param x A list object containing the anomaly detection results from the \code{detectAnomaly} function. 
#' Each element of the list should correspond to a cell type and contain \code{reference_mat_subset}, \code{query_mat_subset}, 
#' \code{var_explained}, and \code{anomaly}.
#' @param cell_type A character string specifying the cell type for which the plots should be generated. This should
#' be a name present in \code{x}. If NULL, the "Combined" cell type will be plotted. Default is NULL.
#' @param pc_subset A numeric vector specifying the indices of the PCs to be included in the plots. If NULL, all PCs
#' in \code{reference_mat_subset} will be included.
#' @param data_type A character string specifying whether to plot the "query" data or the "reference" data. Default is "query".
#' @param ... Additional arguments.
#' 
#' @return A ggplot2 object representing the PCA plots with anomalies highlighted.
#' 
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{detectAnomaly}}
#' 
#' @rdname detectAnomaly
#' 
# Function to create faceted scatter plots for specified PC combinations
plot.detectAnomaly <- function(x, cell_type = NULL, pc_subset = NULL, data_type = c("query", "reference"), ...) {
    
    # Check if PCA was used for computations
    if(!("var_explained" %in% names(x[[names(x)[1]]])))
        stop("The plot function can only be used if \'n_components\' is not NULL.")
    
    # Check input for cell type
    if(is.null(cell_type)){
        cell_type <- "Combined"
    } else{
        if(!(cell_type %in% names(x)))
            stop("\'cell_type\' is not available in \'x\'.")
    }
    
    # Check input for pc_subset
    if(!is.null(pc_subset)){
        if(!all(pc_subset %in% 1:ncol(x[[cell_type]]$reference_mat_subset)))
            stop("\'pc_subset\' is out of range.")
    } else{
        pc_subset <- 1:ncol(x[[cell_type]]$reference_mat_subset)
    }
    
    # Check input for data_type
    data_type <- match.arg(data_type)
    
    # Filter data to include only specified PCs
    if(is.null(x[[cell_type]]$query_mat_subset) && data_type == "query"){
        stop("There is no query data available in the \'detectAnomaly\' object.")
    } else{
        if(data_type == "query"){
            data_subset <- x[[cell_type]]$query_mat_subset[, pc_subset, drop = FALSE]
            anomaly <- x[[cell_type]]$query_anomaly
            
        } else if(data_type == "reference"){
            data_subset <- x[[cell_type]]$reference_mat_subset[, pc_subset, drop = FALSE]
            anomaly <- x[[cell_type]]$reference_anomaly
        }
    }
    
    # Modify column names to include percentage of variance explained
    colnames(data_subset) <- paste0("PC", pc_subset, 
                                    " (", sprintf("%.1f%%", x[[cell_type]]$var_explained[pc_subset] * 100), ")")
    
    # Create all possible pairs of specified PCs
    pc_names <- colnames(data_subset)
    pairs <- expand.grid(x = pc_names, y = pc_names)
    pairs <- pairs[pairs$x != pairs$y, ]
    
    # Create a new data frame with all possible pairs of specified PCs
    data_pairs_list <- lapply(1:nrow(pairs), function(i) {
        x_col <- pairs$x[i]
        y_col <- pairs$y[i]
        data_frame <- data.frame(data_subset[, c(x_col, y_col)])
        colnames(data_frame) <- c("x_value", "y_value")
        data_frame$x <- x_col
        data_frame$y <- y_col
        data_frame
    })
    data_pairs <- do.call(rbind, data_pairs_list)
    
    # Remove redundant data (to avoid duplicated plots)
    data_pairs <- data_pairs[as.numeric(data_pairs$x) < as.numeric(data_pairs$y),]
    
    # Add anomalies vector to data_pairs dataframe
    data_pairs$anomaly <- rep(anomaly, choose(length(pc_subset), 2))
    
    # Create the ggplot object with facets
    plot <- ggplot2::ggplot(data_pairs, ggplot2::aes(x = x_value, y = y_value, color = factor(anomaly))) +
        ggplot2::geom_point(size = 1, alpha = 0.5) + 
        ggplot2::scale_color_manual(values = c("black", "red"), labels = c("Normal", "Anomaly")) + 
        ggplot2::facet_grid(rows = ggplot2::vars(y), cols = ggplot2::vars(x), scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10)) + 
        ggplot2::labs(title = paste0("Isolation Forest Anomaly Plot: ", cell_type), color = "iForest Type")
    print(plot)
}
