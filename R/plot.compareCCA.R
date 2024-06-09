#' @title Plot Visualization of Output from compareCCA Function
#' 
#' @description This function generates a visualization of the output from the `compareCCA` function.
#' The plot shows the cosine similarities of canonical correlation analysis (CCA) coefficients,
#' with point sizes representing the correlations.
#'
#' @details The function converts the input list into a data frame suitable for plotting with `ggplot2`.
#' Each point in the scatter plot represents the cosine similarity of CCA coefficients, with the size of the point
#' indicating the correlation.
#'
#' @param x A list containing the output from the `compareCCA` function. 
#' This list should include `cosine_similarity` and `correlations`.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the scatter plot of cosine similarities of CCA coefficients and correlations.
#'
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{compareCCA}}
#' 
#' @rdname compareCCA
#' 
# Plot visualization of output from compareCCA function
plot.compareCCA <- function(x, ...){
    
    # Create a data frame for plotting
    comparison_data <- data.frame(CCA = paste0("CC", 1:length(x$correlations)),
                                  Cosine = x$cosine_similarity,
                                  Correlation = x$correlations)
    comparison_data$CC <- factor(comparison_data$CCA, levels = comparison_data$CCA)
    
    
    cca_plot <- ggplot2::ggplot(comparison_data, aes(x = CCA, y = Cosine, size = Correlation)) +
        ggplot2::geom_point() +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(title = "Cosine Similarities of CCA Coefficients with Correlation",
                      x = "",
                      y = "Cosine Similarity of CC Coefficients",
                      size = "Correlation") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, 
                                                           size = 12, hjust = 1))
    return(cca_plot)
}