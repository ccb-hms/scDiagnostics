#' @title Plot Wasserstein Distance Distributions for Multiple Cell Types
#'
#' @description
#' This function generates ridge plots comparing reference-reference
#' and reference-query Wasserstein distance distributions for each cell type,
#' along with the probability of superiority metric.
#'
#' @details
#' The function creates faceted ridge plots showing two clearly separated density curves for each
#' cell type: one for the reference-reference distribution (null) and one for the
#' reference-query distribution. The probability of superiority is displayed for each
#' comparison, representing the probability that a ref-query distance exceeds a ref-ref distance.
#'
#' @param x A list object containing the Wasserstein distance results from the \code{calculateWassersteinDistance} function.
#' @param cell_types A character vector specifying which cell types to plot. If NULL, all cell types are plotted.
#' @param ... Additional arguments for future extensions.
#'
#' @keywords internal
#'
#' @return A ggplot2 object representing the comparison of Wasserstein distance distributions.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateWassersteinDistance}}
#'
#' @rdname calculateWassersteinDistance
#'
# Function to plot densities of Wasserstein distances
plot.calculateWassersteinDistanceObject <- function(
        x,
        cell_types = NULL,
        ...){

    # Determine which cell types to plot
    if (is.null(cell_types)) {
        cell_types <- x$cell_types
    } else {
        # Check that requested cell types are available
        missing_types <- setdiff(cell_types, x$cell_types)
        if (length(missing_types) > 0) {
            warning(paste("Cell types not found in data:",
                          paste(missing_types, collapse = ", ")))
        }
        cell_types <- intersect(cell_types, x$cell_types)
    }

    if (length(cell_types) == 0) {
        stop("No valid cell types to plot.")
    }

    # Prepare data for plotting
    plot_data <- data.frame()

    for (cell_type in cell_types) {
        # Reference-reference distribution
        ref_ref_data <- data.frame(
            wasserstein_dist = x$ref_ref_dist[[cell_type]],
            distribution = "Reference-Reference",
            cell_type = cell_type
        )

        # Reference-query distribution
        ref_query_data <- data.frame(
            wasserstein_dist = x$ref_query_dist[[cell_type]],
            distribution = "Reference-Query",
            cell_type = cell_type
        )

        plot_data <- rbind(plot_data, ref_ref_data, ref_query_data)
    }

    # Set factor levels to control order (Reference-Reference at bottom, Reference-Query at top)
    plot_data$distribution <- factor(plot_data$distribution,
                                     levels = c("Reference-Reference",
                                                "Reference-Query"))

    # Set factor levels for cell_type to preserve the order specified in cell_types parameter
    plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_types)

    # Setting up color data
    plot_data[["cell_type_distribution"]] <-
        paste(plot_data[["cell_type"]],
              plot_data[["distribution"]])
    plot_data[["cell_type_distribution"]] <-
        factor(plot_data[["cell_type_distribution"]],
               levels = unique(plot_data[["cell_type_distribution"]]))
    cell_type_distribution_colors <-
        generateColors(levels(plot_data[["cell_type_distribution"]]),
                       paired = TRUE)

    # Prepare annotation data
    annotation_data <- data.frame(
        cell_type = names(x$probability_superiority),
        probability = x$probability_superiority,
        label = paste0("P(Ref-Query > Ref-Ref): ",
                       sprintf("%.3f", x$probability_superiority))
    )

    # Filter annotation data to match selected cell types and preserve order
    annotation_data <- annotation_data[annotation_data$cell_type %in% cell_types, ]
    annotation_data$cell_type <- factor(annotation_data$cell_type, levels = cell_types)

    # Create the ridge plot
    ridge_plot <- ggplot2::ggplot(plot_data,
                                  ggplot2::aes(
                                      x = .data[["wasserstein_dist"]],
                                      y = .data[["distribution"]],
                                      fill = .data[["cell_type_distribution"]])) +
        ggridges::geom_density_ridges(
            alpha = 0.7,
            scale = 0.8,
            rel_min_height = 0.01,
            jittered_points = FALSE,
            position = "identity"
        ) +
        ggplot2::facet_wrap(~ .data[["cell_type"]], scales = "free_x", ncol = 2) +
        ggplot2::scale_fill_manual(
            name = "Distribution",
            values = cell_type_distribution_colors,
            labels = names(cell_type_distribution_colors)
        ) +
        ggplot2::scale_y_discrete(
            expand = ggplot2::expansion(mult = c(0.02, 0.02))
        ) +
        ggplot2::labs(
            title = "Comparison of Wasserstein Distance Distributions by Cell Type",
            x = "Wasserstein Distance",
            y = "Distribution"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.border = ggplot2::element_rect(color = "black",
                                                 fill = NA,
                                                 linewidth = 0.5),
            axis.title.x = ggplot2::element_text(size = 12),
            legend.position = "bottom",
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.ticks.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold", hjust = 0.5),
            strip.text = ggplot2::element_text(size = 11,
                                               face = "bold"),
            strip.background = ggplot2::element_rect(fill = "white",
                                                     color = "black",
                                                     linewidth = 0.5),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(color = "gray",
                                                       linetype = "dotted"),
            panel.grid.major.y = ggplot2::element_blank()
        ) +
        ggplot2::geom_text(
            data = annotation_data,
            ggplot2::aes(
                x = Inf, y = Inf,
                label = .data[["label"]]
            ),
            hjust = 1.05, vjust = 1.1,
            inherit.aes = FALSE,
            size = 3,
            fontface = "bold",
            color = "darkblue"
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(override.aes = list(alpha = 0.8))
        )

    return(ridge_plot)
}
