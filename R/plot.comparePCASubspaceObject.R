#' @title Plot Visualization of Output from comparePCASubspace Function
#'
#' @description
#' The S3 plot method generates a visualization of the output from the \code{comparePCASubspace} function.
#' The plot shows the cosine of principal angles between reference and query principal components,
#' with point sizes representing the variance explained and colors showing the difference in variance between datasets.
#'
#' @details
#' The S3 plot method converts the input list into a data frame suitable for plotting with \code{ggplot2}.
#' Each point in the scatter plot represents the cosine of a principal angle, with the size of the point
#' indicating the average variance explained by the corresponding principal components. The color represents
#' the difference in variance explained between reference and query datasets.
#'
#' @param x A numeric matrix output from the \code{comparePCASubspace} function, representing
#' cosine similarities between query and reference principal components.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @keywords internal
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the cosine similarities with variance information.
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
plot.comparePCASubspaceObject <- function(x,
                                          ...){

    # Create a data frame for plotting
    plot_data <- data.frame(
        PC = paste0("Ref PC", x[["cosine_id"]][, 1], " - Query PC", x[["cosine_id"]][, 2]),
        PC_Label = paste0("Ref PC", x[["cosine_id"]][, 1],
                          " (", round(x[["var_explained_ref"]], 1), "%) - Query PC",
                          x[["cosine_id"]][, 2], " (", round(x[["var_explained_query"]], 1), "%)"),
        Cosine = x[["cosine_similarity"]],
        VarianceExplained = x[["var_explained_avg"]],
        RefVariance = x[["var_explained_ref"]],
        QueryVariance = x[["var_explained_query"]],
        VarianceDiff = x[["var_explained_ref"]] - x[["var_explained_query"]]
    )
    plot_data[["PC_Label"]] <- factor(plot_data[["PC_Label"]], levels = plot_data[["PC_Label"]])

    # Create plot
    pc_plot <- ggplot2::ggplot(
        plot_data, ggplot2::aes(x = .data[["PC_Label"]],
                                y = .data[["Cosine"]],
                                size = .data[["VarianceExplained"]],
                                color = .data[["VarianceDiff"]])) +
        ggplot2::geom_point(alpha = 0.8, stroke = 1) +
        ggplot2::scale_size_continuous(range = c(3, 10), name = "Avg Variance\nExplained (%)") +
        ggplot2::scale_color_gradient2(
            low = "blue", mid = "gray", high = "red", midpoint = 0,
            name = "Variance Diff\n(Ref - Query)"
        ) +
        ggplot2::labs(
            title = "Principal Angles Cosines with Variance Explained",
            subtitle = paste0("Weighted Cosine Similarity: ", sprintf("%.3f", x[["weighted_cosine_similarity"]])),
            x = "",
            y = "Cosine Similarity of Principal Angle"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                size = 10, hjust = 1),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5)
        )
    return(pc_plot)
}
