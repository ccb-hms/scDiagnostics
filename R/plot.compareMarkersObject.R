#' @title Plot Visualization of Output from `compareMarkers` Function
#'
#' @description
#' The S3 plot method generates a comprehensive visualization of the output from the `compareMarkers` function.
#' The plot shows marker gene overlap and expression consistency between query and reference cell types,
#' with quality assessment and detailed annotations.
#'
#' @details
#' The S3 plot method creates a scatter plot showing the relationship between marker overlap (x-axis)
#' and expression consistency (y-axis) for each cell type. Points are colored by quality score and
#' sized by the minimum number of cells. Quality zones provide visual guidance for interpretation.
#'
#' @param x A list containing the output from the \code{compareMarkers} function.
#' @param cell_types Character vector specifying which cell types to plot. If NULL, all cell types are plotted.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @keywords internal
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the marker gene comparison results.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{compareMarkers}}
#'
#' @rdname compareMarkers
#'
# Plot visualization of output from compareMarkers function
plot.compareMarkersObject <- function(x, cell_types = NULL, ...){

    # Filter data by cell types if specified
    if (!is.null(cell_types)) {
        available_types <- intersect(cell_types, names(x[["marker_overlap"]]))
        if (length(available_types) == 0) {
            stop("None of the specified cell types are available in the data")
        }

        # Filter all relevant components
        x[["marker_overlap"]] <-
            x[["marker_overlap"]][available_types]
        x[["expression_consistency"]] <-
            x[["expression_consistency"]][available_types]
        x[["quality_scores"]] <-
            x[["quality_scores"]][available_types]
        x[["n_cells_query"]] <-
            x[["n_cells_query"]][available_types]
        x[["n_cells_ref"]] <-
            x[["n_cells_ref"]][available_types]
        x[["selected_cell_types"]] <- available_types
    }

    # Create data frame for plotting
    plot_data <- data.frame(
        Cell_Type = as.character(names(x[["marker_overlap"]])),
        Marker_Overlap = as.numeric(x[["marker_overlap"]]),
        Expression_Consistency = as.numeric(x[["expression_consistency"]]),
        Quality = as.character(x[["quality_scores"]]),
        stringsAsFactors = FALSE
    )

    # Handle cell counts properly
    query_counts <- x[["n_cells_query"]]
    ref_counts <- x[["n_cells_ref"]]

    # Ensure we have counts for all cell types
    plot_data[["N_Query"]] <-
        as.numeric(query_counts[plot_data[["Cell_Type"]]])
    plot_data[["N_Ref"]] <-
        as.numeric(ref_counts[plot_data[["Cell_Type"]]])

    # Calculate minimum cell count
    plot_data[["Min_Cells"]] <- pmin(plot_data[["N_Query"]],
                                     plot_data[["N_Ref"]])

    # Define colors
    quality_colors <- c("Good" = "#2E8B57",
                        "Moderate" = "#FF8C00",
                        "Poor" = "#DC143C")

    # Create formal filter description
    filter_description <- switch(
        x[["anomaly_filter_used"]],
        "none" = "",
        "anomalous_only" = " (Anomalous Cells Only)",
        "non_anomalous_only" = " (Non-Anomalous Cells Only)")

    # Create subtitle with filtering information
    subtitle_text <- paste0("Analysis of ",
                            length(x[["selected_cell_types"]]),
                            " cell types",
                            filter_description)

    # Create the main scatter plot
    marker_plot <- ggplot2::ggplot(plot_data,
                                   ggplot2::aes(
                                       x = .data[["Marker_Overlap"]],
                                       y = .data[["Expression_Consistency"]])) +

        # Add quality zone backgrounds
        ggplot2::annotate("rect", xmin = 0.7, xmax = 1,
                          ymin = 0.7, ymax = 1,
                          fill = "#2E8B57", alpha = 0.1) +
        ggplot2::annotate("rect", xmin = 0.4, xmax = 0.7,
                          ymin = 0.4, ymax = 0.7,
                          fill = "#FF8C00", alpha = 0.1) +
        ggplot2::annotate("rect", xmin = 0, xmax = 0.4,
                          ymin = 0, ymax = 0.4,
                          fill = "#DC143C", alpha = 0.1) +

        # Add threshold lines
        ggplot2::geom_hline(yintercept = 0.7,
                            linetype = "dashed",
                            color = "#2E8B57",
                            alpha = 0.5) +
        ggplot2::geom_vline(xintercept = 0.7,
                            linetype = "dashed",
                            color = "#2E8B57",
                            alpha = 0.5) +
        ggplot2::geom_hline(yintercept = 0.4,
                            linetype = "dashed",
                            color = "#FF8C00",
                            alpha = 0.5) +
        ggplot2::geom_vline(xintercept = 0.4,
                            linetype = "dashed",
                            color = "#FF8C00",
                            alpha = 0.5) +

        # Add points
        ggplot2::geom_point(ggplot2::aes(
            color = .data[["Quality"]],
            size = .data[["Min_Cells"]]),
            alpha = 0.8) +

        # Add text labels for cell types
        ggplot2::geom_text(ggplot2::aes(
            label = .data[["Cell_Type"]]),
            size = 3,
            hjust = -0.1,
            vjust = -0.1,
            alpha = 0.8) +

        # Color and size scales
        ggplot2::scale_color_manual(values = quality_colors,
                                    name = "Quality\n(Worst of Both)") +
        ggplot2::scale_size_continuous(range = c(3, 8),
                                       name = "Min Cells\nper Type") +

        # Set axis limits
        ggplot2::xlim(0, 1) +
        ggplot2::ylim(0, 1) +

        # Labels and theme
        ggplot2::labs(
            title = "Marker Gene Comparison Between Query and Reference",
            subtitle = subtitle_text,
            x = "Marker Gene Overlap (Jaccard Index)",
            y = "Expression Consistency of Reference Markers"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            legend.position = "right"
        )

    # Add text annotations for zones
    marker_plot <- marker_plot +
        ggplot2::annotate("text", x = 0.85, y = 0.95, label = "Good",
                          size = 4, color = "#2E8B57", fontface = "bold") +
        ggplot2::annotate("text", x = 0.55, y = 0.55, label = "Moderate",
                          size = 4, color = "#FF8C00", fontface = "bold") +
        ggplot2::annotate("text", x = 0.2, y = 0.2, label = "Poor",
                          size = 4, color = "#DC143C", fontface = "bold")

    return(marker_plot)
}
