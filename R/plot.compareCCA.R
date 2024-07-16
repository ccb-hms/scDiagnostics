#' @title Plot Visualization of Output from compareCCA Function
#' 
#' @description 
#' The S3 plot method generates a visualization of the output from the `compareCCA` function.
#' The plot shows the cosine similarities of canonical correlation analysis (CCA) coefficients,
#' with point sizes representing the correlations.
#'
#' @details 
#' The S3 plot method converts the input list into a data frame suitable for plotting with \code{ggplot2}.
#' Each point in the scatter plot represents the cosine similarity of CCA coefficients, with the size of the point
#' indicating the correlation.
#'
#' @param x A list containing the output from the \code{compareCCA} function. 
#' This list should include \code{cosine_similarity} and \code{correlations}.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the scatter plot of cosine similarities of CCA coefficients and correlations.
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
    comparison_data <- data.frame(CCA = paste0("CC", seq_len(length(x[["correlations"]]))),
                                  Cosine = x[["cosine_similarity"]],
                                  Correlation = x[["correlations"]])
    comparison_data$CC <- factor(comparison_data[["CCA"]], levels = comparison_data[["CCA"]])
    
    
    cca_plot <- ggplot2::ggplot(comparison_data, ggplot2::aes(x = .data[["CCA"]], y = .data[["Cosine"]], 
                                                              size = .data[["Correlation"]])) +
        ggplot2::geom_point() +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(title = "Cosine Similarities of CCA Coefficients with Correlation",
                      x = "", y = "Cosine Similarity of CC Coefficients", size = "Correlation") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
                       axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"))
    return(cca_plot)
}
