#' @title Plot Density of Probabilities for Cell Type Classification
#'
#' @description 
#' The S3 plot method generates a density plot showing the distribution of probabilities for each cell of belonging to 
#' either the reference or query dataset for each cell type.
#'
#' @details 
#' The S3 plot method creates a density plot to visualize the distribution of probabilities for each cell belonging to the 
#' reference or query dataset for each cell type. It utilizes the ggplot2 package for plotting.
#'
#' @param x An object of class \code{nearestNeighbotDiagnostics} containing the probabilities calculated by the \code{\link{calculateNearestNeighborProbabilities}} function.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types in \code{x} will be plotted. Default is NULL.
#' @param ... Additional arguments to be passed to \code{\link[ggplot2]{geom_density}}.
#' 
#' @keywords internal
#'
#' @return The S3 plot method returns a \code{ggplot} density plot.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{calculateNearestNeighborProbabilities}}
#' 
#' @rdname calculateNearestNeighborProbabilities
#' 
# Function to plot probabilities of each cell of belonging to reference or query dataset for each cell type
plot.calculateNearestNeighborProbabilities <- function(x, cell_types = NULL, ...) {
    
    if(!is.null(cell_types)){
        
        if(!all(cell_types %in% names(x)))
            stop("One or more of the \'cell_types'\ is not available.")
        
        # Subset cell types
        x <- x[cell_types]
    }

    # Create a data frame with the x values, densities, and quantiles
    densities <- do.call(rbind, lapply(seq_len(length(x)), function(i) {
        x_range <- seq(0.5 - 3 * sqrt(0.25 / x[[i]][["n_neighbor"]] / x[[i]][["n_query"]]), 
                       0.5 + 3 * sqrt(0.25 / x[[i]][["n_neighbor"]] / x[[i]][["n_query"]]), length.out = 1000)
        data.frame(x = x_range, y = dnorm(x_range, mean = 0.5, sd = sqrt(0.25 / x[[i]][["n_neighbor"]] / x[[i]][["n_query"]])), 
                   cell_type = names(x)[i],
                   query_prob = x[[i]][["query_prob"]],
                   lb = qnorm(0.025, mean = 0.5, sd = sqrt(0.25 / x[[i]][["n_neighbor"]] / x[[i]][["n_query"]])),
                   ub = qnorm(0.975, mean = 0.5, sd = sqrt(0.25 / x[[i]][["n_neighbor"]] / x[[i]][["n_query"]])))
    }))
    
    # Define the colors for cell types
    cell_type_colors <- generateColors(sort(unique(densities[["cell_type"]])), paired = FALSE)
    
    # Plot the density functions with facets and quantile lines
    density_plot <- ggplot2::ggplot(densities, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], 
                                                            fill = .data[["cell_type"]])) +
        ggplot2::geom_line(color = "black") +
        ggplot2::geom_area(alpha = 0.7) + 
        ggplot2::facet_wrap(~ .data[["cell_type"]], scales = "free") +
        ggplot2::scale_fill_manual(values = cell_type_colors, name = "Cell Type") + 
        ggplot2::geom_vline(data = densities, ggplot2::aes(xintercept = .data[["query_prob"]]), linetype = "dashed", color = "black") +
        ggplot2::geom_vline(data = densities, ggplot2::aes(xintercept = .data[["lb"]]), linetype = "solid", color = "black") +
        ggplot2::geom_vline(data = densities, ggplot2::aes(xintercept = .data[["ub"]]), linetype = "solid", color = "black") +
        ggplot2::labs(title = "Approximate Density Functions Under Null Hypothesis (Equal Reference and Query Distributions)", 
                      x = "Probability Query Cell", y = "Density") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10))
    return(density_plot)
}
