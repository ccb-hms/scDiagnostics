#' @title Plot Visualization of Output from comparePCASubspace Function
#' 
#' @description This function generates a visualization of the output from the `comparePCASubspace` function.
#' The plot shows the cosine of principal angles between reference and query principal components,
#' with point sizes representing the variance explained.
#' 
#' @details The function converts the input list into a data frame suitable for plotting with `ggplot2`.
#' Each point in the scatter plot represents the cosine of a principal angle, with the size of the point
#' indicating the average variance explained by the corresponding principal components.
#' 
#' @param x A numeric matrix output from the `comparePCA` function, representing 
#' cosine similarities between query and reference principal components.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the heatmap of cosine similarities.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{comparePCASubspace}}
#' 
#' @rdname comparePCASubspace
#'
# Function to produce the visualization of output from comparePCASubspace function
plot.comparePCASubspace <- function(x, ...){
    
    # Create a data frame for plotting
    x <- data.frame(PC = paste0("Ref PC", x$cosine_id[, 1],
                                " - Query PC", x$cosine_id[, 2]),
                    Cosine = x$cosine_similarity,
                    VarianceExplained = subspace_comparison$var_explained_avg)
    x$PC <- factor(x$PC, levels = x$PC)
    
    # Create plot
    pc_plot <- ggplot2::ggplot(x, aes(x = PC, y = Cosine, size = VarianceExplained)) +
        ggplot2::geom_point() +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(title = "Principal Angles Cosines with Variance Explained",
                      x = "",
                      y = "Cosine Similarity of Principal Angle",
                      size = "Variance Explained") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, 
                                                           size = 12, hjust = 1))
    return(pc_plot)
}