#' @title 
#' Plot the output of the calculateAveragePairwiseCorrelation function
#'
#' @description 
#' This function takes the output of the calculateAveragePairwiseCorrelation function,
#' which should be a matrix of pairwise correlations, and plots it as a heatmap.
#' 
#' @details 
#' This function converts the correlation matrix into a dataframe, creates a heatmap using ggplot2,
#' and customizes the appearance of the heatmap with updated colors and improved aesthetics.
#'
#' @param x Output matrix from calculateAveragePairwiseCorrelation function.
#' @param ... Additional arguments to be passed to the plotting function.
#'
#' @return A ggplot2 object representing the heatmap plot.
#' 
#' @export
#'         
#' @seealso \code{\link{calculateAveragePairwiseCorrelation}}
#' 
#' @rdname calculateAveragePairwiseCorrelation
#' 
# Function to plot the output of the calculateAveragePairwiseCorrelation function
plot.calculateAveragePairwiseCorrelation <- function(x, ...){
    
    # Convert matrix to dataframe
    cor_df <- as.data.frame(as.table(x))
    cor_df$Var1 <- factor(cor_df$Var1, levels = rownames(x))
    cor_df$Var2 <- factor(cor_df$Var2, levels = rev(colnames(x)))
    
    # Create the heatmap with updated colors and improved aesthetics
    heatmap_plot <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var2, y = Var1)) +
        ggplot2::geom_tile(ggplot2::aes(fill = Freq), color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = round(Freq, 2)), color = "black", size = 3, family = "sans") +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                                      midpoint = 0, limits = c(min(cor_df$Freq), max(cor_df$Freq)),
                                      name = "Correlation",
                                      breaks = seq(-1, 1, by = 0.2)) +  # Specify color scale breaks
        ggplot2::labs(title = "Correlation Heatmap", x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
                       axis.text.y = ggplot2::element_text(family = "sans"),  # Set font family for y-axis labels
                       plot.title = ggplot2::element_text(face = "bold"),  # Make title bold
                       legend.position = "right",  # Place legend on RHS
                       legend.title = ggplot2::element_text(face = "italic"))
    
    # Print the plot
    print(heatmap_plot)
}
