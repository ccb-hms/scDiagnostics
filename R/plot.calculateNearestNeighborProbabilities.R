#' @title Plot Density of Probabilities for Cell Type Classification
#'
#' @description This function generates a density plot showing the distribution of probabilities for each sample of belonging to 
#' either the reference or query dataset for each cell type.
#'
#' @details This function creates a density plot to visualize the distribution of probabilities for each sample belonging to the 
#' reference or query dataset for each cell type. It utilizes the ggplot2 package for plotting.
#'
#' @param x An object of class \code{nearestNeighbotDiagnostics} containing the probabilities calculated by the \code{\link{calculateNearestNeighborProbabilities}} function.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types in \code{x} will be plotted. Default is NULL.
#' @param prob_type A character string specifying the type of probability to plot. Must be either "query" or "reference". Default is "query".
#' @param ... Additional arguments to be passed to \code{\link[ggplot2]{geom_density}}.
#'
#' @return A ggplot2 density plot.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{calculateNearestNeighborProbabilities}}
#' 
#' @rdname calculateNearestNeighborProbabilities
#' 
# Function to plot probabilities of each sample of belonging to reference or query dataset for each cell type
plot.calculateNearestNeighborProbabilities <- function(x, cell_types = NULL,
                                                       prob_type = c("query", "reference")[1], ...) {
    
    # Check input for probability type
    if(!(prob_type %in% c("query", "reference")))
        stop("\'prob_type\' must be one of \'query\' or \'reference\'.")
    
    # Convert probabilities to data frame
    probabilities_df <- do.call(rbind, lapply(names(x), function(ct) {
        data.frame(cell_types = ct, 
                   probability = x[[ct]][[ifelse(prob_type == "reference", "prob_ref", "prob_query")]])
    }))
    
    if(!is.null(cell_types)){
        
        if(!all(cell_types %in% names(x)))
            stop("One or more of the \'cell_types'\ is not available.")
        
        # Subset cell types
        probabilities_df <- probabilities_df[probabilities_df$cell_types %in% cell_types,]
    }

    # Create density plot
    density_plot <- ggplot2::ggplot(probabilities_df, ggplot2::aes(x = probability, fill = cell_types)) +
        ggplot2::geom_density(alpha = 0.7) +
        ggplot2::labs(x = "Probability", y = "Density", title = "Density Plot of Probabilities") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10)) +
        ggplot2::facet_wrap(~cell_types, scales = "free")
    if(length(unique(probabilities_df[["cell_types"]])) > 1){
        
        # Setting up colors
        paired_palette <- RColorBrewer::brewer.pal(12, "Paired")
        dark_palette <- paired_palette[seq(2, length(paired_palette), by = 2)]
        order_combinations <- sort(unique(probabilities_df[["cell_types"]]))
        dark_colors_named <- setNames(dark_palette[1:length(order_combinations)], order_combinations)
        density_plot <- density_plot + 
            ggplot2::scale_fill_manual(values = dark_colors_named)
    }
    
    return(density_plot)
}
