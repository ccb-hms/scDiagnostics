#' @title Plot Projected Data on Discriminant Spaces
#'
#' @description
#' The S3 plot method visualizes the projected reference and query data on discriminant spaces using either a scatterplot, boxplot, or varplot.
#'
#' @details
#' - **Scatterplot**: Displays projected data points, with colors used to differentiate between cell types and datasets.
#' - **Boxplot**: Shows the distribution of projected data values for each cell type, separated by datasets.
#' - **Varplot**: Highlights the top contributing variables for each discriminant axis, differentiating between positive and negative loadings.
#'
#' @param x An object of class \code{calculateSIRSpace}, which contains projected data on the discriminant space.
#' Each element should include 'ref_proj' and 'query_proj' data frames representing reference and query projections.
#' @param plot_type A character string indicating the type of plot to generate. Options are "scatterplot", "boxplot", or "varplot". Default is "scatterplot".
#' @param sir_subset A numeric vector specifying which discriminant axes (SIR components) to include in the plot. Default is the first 5 axes.
#' @param n_top_vars Number of top contributing variables to display in varplot. Default is 5.
#' @param ... Additional arguments to be passed to the plotting functions.
#'
#' @keywords internal
#'
#' @return A \code{ggplot} object representing the chosen visualization (scatterplot, boxplot, or varplot) of the projected data.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateSIRSpace}}
#'
#' @rdname calculateSIRSpace
#'
#' @importFrom utils head
#'
# This function plots the SIR output either as a scatterplot for a boxplot
plot.calculateSIRSpace <- function(x,
                                   plot_type = c("scatterplot", "boxplot", "varplot"),
                                   sir_subset = NULL,
                                   n_top_vars = 5,
                                   ...){

    # Match argument for plot type
    plot_type <- match.arg(plot_type)

    # Value for SIR subset
    if(is.null(sir_subset)){
        sir_subset <- seq_len(min(5, ncol(x[["rotation_mat"]])))
    }

    # Check input for SIR subset
    if(any(!(sir_subset %in% seq_len(ncol(x[["rotation_mat"]]))))){
        stop("\"sir_subset\" is out of range.")
    }

    # Identify cell types
    cell_types <- unique(rownames(x[["cond_means"]]))

    if(plot_type == "scatterplot"){

        # Remove cells with NA cell type
        sir_projections <- na.omit(x[["sir_projections"]])

        # Create a subset of plot names for the SIR subset
        plot_names <- paste0(
            "SIR", sir_subset, " (",
            sprintf("%.1f%%", x[["percent_var"]][sir_subset]), ")")

        # Create all possible pairs of the specified PCs in the subset
        pairs <- expand.grid(x = sir_subset, y = sir_subset)
        pairs <- pairs[pairs[["x"]] < pairs[["y"]], ]

        # Create a new data frame with all possible pairs of specified PCs
        data_pairs_list <- lapply(seq_len(nrow(pairs)), function(i) {
            x_col <- pairs[["x"]][i]
            y_col <- pairs[["y"]][i]

            # Properly subset the projection data using sir_subset indices
            data_frame <- data.frame(
                sir_projections[, c(x_col, y_col)],
                paste(sir_projections[["dataset"]], sir_projections[["cell_type"]], sep = " "))

            colnames(data_frame) <- c("x_value", "y_value", "cell_type_dataset")
            data_frame[["x"]] <- plot_names[match(x_col, sir_subset)]
            data_frame[["y"]] <- plot_names[match(y_col, sir_subset)]

            return(data_frame)
        })

        # Combine all pairs into one data frame
        data_pairs <- do.call(rbind, data_pairs_list)
        data_pairs[["x"]] <- factor(data_pairs[["x"]],
                                    levels = plot_names)
        data_pairs[["y"]] <- factor(data_pairs[["y"]],
                                    levels = plot_names)

        # Define the order of cell type and dataset combinations
        order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)), rep(sort(cell_types), each = 2))
        data_pairs[["cell_type_dataset"]] <- factor(data_pairs[["cell_type_dataset"]], levels = order_combinations)
        cell_type_colors <- generateColors(order_combinations, paired = TRUE)

        # Set SIR identity as factor
        data_pairs[["x"]] <- factor(data_pairs[["x"]], levels = unique(data_pairs[["x"]]))
        data_pairs[["y"]] <- factor(data_pairs[["y"]], levels = unique(data_pairs[["y"]]))

        # Create the ggplot object (with facets if PCA)
        plot_obj <- ggplot2::ggplot(
            data_pairs, ggplot2::aes(x = .data[["x_value"]],
                                     y = .data[["y_value"]],
                                     color = .data[["cell_type_dataset"]])) +
            ggplot2::geom_point(alpha = 0.5, size = 1) +
            ggplot2::scale_color_manual(values = cell_type_colors,
                                        name = "Cell Types") +
            ggplot2::facet_grid(rows = ggplot2::vars(.data[["y"]]),
                                cols = ggplot2::vars(.data[["x"]]),
                                scales = "free") +
            ggplot2::xlab("") + ggplot2::ylab("") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_line(color = "gray",
                                                         linetype = "dotted"),
                plot.title = ggplot2::element_text(size = 14,
                                                   face = "bold",
                                                   hjust = 0.5),
                axis.title = ggplot2::element_text(size = 12),
                axis.text = ggplot2::element_text(size = 10))

    } else if (plot_type == "boxplot"){

        # Remove cells with NA cell type
        sir_projections <- na.omit(x[["sir_projections"]])

        # Create the long format data frame manually
        sir_projections <- sir_projections[!is.na(sir_projections[["cell_type"]]),]
        if(!is.null(cell_types)){
            if(all(cell_types %in% sir_projections[["cell_type"]])){
                sir_projections <- sir_projections[which(sir_projections[["cell_type"]] %in%
                                                             cell_types),]
            } else{
                stop("One or more of the specified \'cell_types\' are not available.")
            }
        }
        sir_long <- data.frame(SIR = rep(paste0("sir", sir_subset),
                                         each = nrow(sir_projections)),
                               Value = unlist(c(sir_projections[, sir_subset])),
                               dataset = rep(sir_projections[["dataset"]],
                                             length(sir_subset)),
                               cell_type = rep(sir_projections[["cell_type"]],
                                               length(sir_subset)))
        sir_long[["SIR"]] <- toupper(sir_long[["SIR"]])

        # Create a new variable representing the combination of cell type and dataset
        sir_long[["cell_type_dataset"]] <- paste(sir_long[["dataset"]],
                                                 sir_long[["cell_type"]],
                                                 sep = " ")

        # Define the order of cell type and dataset combinations
        order_combinations <- paste(rep(c("Reference", "Query"),
                                        length(unique(sir_long[["cell_type"]]))),
                                    rep(sort(unique(sir_long[["cell_type"]])),
                                        each = 2))

        # Reorder the levels of cell type and dataset factor
        sir_long[["cell_type_dataset"]] <- factor(sir_long[["cell_type_dataset"]],
                                                  levels = order_combinations)

        # Set SIR identity as factor
        sir_long[["SIR"]] <- factor(sir_long[["SIR"]], levels = paste0("SIR", sir_subset))

        # Define the colors for cell types
        cell_type_colors <- generateColors(order_combinations,
                                           paired = TRUE)

        # Create the ggplot
        plot_obj <- ggplot2::ggplot(sir_long, ggplot2::aes(
            x = .data[["cell_type"]],
            y = .data[["Value"]],
            fill = .data[["cell_type_dataset"]])) +
            ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) +
            ggplot2::facet_wrap(~ .data[["SIR"]], scales = "free") +
            ggplot2::scale_fill_manual(values = cell_type_colors,
                                       name = "Cell Types") +
            ggplot2::labs(x = "", y = "SIR Score") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10))

    } else if (plot_type == "varplot"){

        # Rotation matrix
        rotation_mat <- x[["rotation_mat"]][, sir_subset, drop = FALSE]

        # Create an empty list to store individual ggplot objects
        plot_list <- list()

        # For each SIR vector (each column of rotation_mat)
        for (i in seq_len(ncol(rotation_mat))) {

            sir_vector <- rotation_mat[, i]

            # Ensure that the vector has names
            if (is.null(names(sir_vector))) {
                names(sir_vector) <- rownames(rotation_mat)
            }

            # Find the top 5 positive and top 5 negative markers
            sorted_positive <- sort(sir_vector[sir_vector > 0], decreasing = TRUE)
            sorted_negative <- sort(sir_vector[sir_vector < 0], decreasing = FALSE)

            top_positive <- head(sorted_positive, n_top_vars)
            top_negative <- head(sorted_negative, n_top_vars)

            # Create a data frame for this SIR vector with marker names and loadings
            df <- data.frame(
                marker = c(names(top_positive), names(top_negative)),
                loading = c(top_positive, top_negative),
                SIR = rep(paste0("SIR", i), length(c(top_positive, top_negative))),
                Direction = c(rep("Positive", length(top_positive)), rep("Negative", length(top_negative)))
            )

            # Ensure that 'marker' is treated as a factor, preserving the original order
            df[["marker"]] <- factor(df[["marker"]], levels = df[["marker"]])

            # Set SIR identity as factor
            df[["SIR"]] <- factor(df[["SIR"]], levels = paste0("SIR", sir_subset))

            # Create a plot for this SIR vector
            plot_list[[i]] <- ggplot2::ggplot(df, ggplot2::aes(
                x = .data[["marker"]],
                y = .data[["loading"]],
                fill = .data[["Direction"]])) +
                ggplot2::geom_bar(stat = "identity") +
                ggplot2::facet_wrap(~ .data[["SIR"]], scales = "free_y") +
                ggplot2::coord_flip() +
                ggplot2::labs(x = NULL, y = "SIR Loading", title = NULL) +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    legend.position = "none",
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                    axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                    axis.text.x = ggplot2::element_text(hjust = 1, size = 10))
        }

        # Combine all plots into a grid
        plot_obj <- patchwork::wrap_plots(plot_list, ncol = 3)
    }

    # Return plot object
    return(plot_obj)
}
