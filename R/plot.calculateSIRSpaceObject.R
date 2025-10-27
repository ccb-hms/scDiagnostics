#' @title Plot SIR Components for Different Cell Types
#'
#' @description
#' This function plots the Sliced Inverse Regression (SIR) components for different cell types in query and reference datasets.
#'
#' @details
#' This function visualizes the SIR projections for specified cell types, providing a pairs plot of the SIR components.
#' It offers various visualization options for different facets of the plot including scatter plots, contours, ellipses, and density plots.
#' When plot_type is "loadings", it creates horizontal bar plots showing the n_top contributing variables for each SIR component.
#'
#' @param x An object of class \code{calculateSIRSpaceObject} containing SIR projections.
#' @param plot_type A character string specifying the type of plot. Either "scores" (default) for SIR projections or "loadings" for variable loadings.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included. Only used when plot_type = "scores".
#' @param sir_subset A numeric vector specifying which SIR components to include in the plot. Default is 1:5.
#' @param lower_facet Type of plot to use for the lower panels. Either "scatter" (default), "contour", "ellipse", or "blank". Only used when plot_type = "scores".
#' @param diagonal_facet Type of plot to use for the diagonal panels. Either "ridge" (default), "density", "boxplot" or "blank". Only used when plot_type = "scores".
#' @param upper_facet Type of plot to use for the upper panels. Either "blank" (default), "scatter", "contour", or "ellipse". Only used when plot_type = "scores".
#' @param n_top A numeric value specifying the number of n_top variables (by absolute loading value) to display. Default is 10 Only used when plot_type = "loadings".
#' @param max_cells_ref Maximum number of reference cells to include in the plot. If NULL,
#' all available reference cells are plotted. Default is NULL. Only used when plot_type = "scores".
#' @param max_cells_query Maximum number of query cells to include in the plot. If NULL,
#' all available query cells are plotted. Default is NULL. Only used when plot_type = "scores".
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggmatrix object representing a pairs plot of specified SIR components for the given cell types and datasets when plot_type = "scores", or a ggplot object showing loadings when plot_type = "loadings".
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateSIRSpace}}
#'
#' @rdname calculateSIRSpace
#'
# Function to plot data projected onto SIR space
plot.calculateSIRSpaceObject <- function(x,
                                         plot_type = c("scores", "loadings"),
                                         cell_types = NULL,
                                         sir_subset = 1:5,
                                         lower_facet = c("scatter", "contour", "ellipse", "blank"),
                                         diagonal_facet = c("ridge", "density", "boxplot", "blank"),
                                         upper_facet = c("blank", "scatter", "contour", "ellipse"),
                                         n_top = 10,
                                         max_cells_ref = NULL,
                                         max_cells_query = NULL,
                                         ...) {

    # Match arguments
    plot_type <- match.arg(plot_type)
    lower_facet <- match.arg(lower_facet)
    diagonal_facet <- match.arg(diagonal_facet)
    upper_facet <- match.arg(upper_facet)

    # Check sir_subset against available SIR components
    if (any(!(sir_subset %in% seq_len(ncol(x[["rotation_mat"]]))))) {
        stop("sir_subset contains values outside the range of available SIR components")
    }

    # Check n_top parameter for loadings plot
    if (plot_type == "loadings") {
        if (!is.numeric(n_top) || n_top <= 0 || n_top != as.integer(n_top)) {
            stop("n_top must be a positive integer.")
        }
    }

    # Validate max_cells parameters
    if (!is.null(max_cells_ref)) {
        if (!is.numeric(max_cells_ref) || max_cells_ref <= 0 || max_cells_ref != as.integer(max_cells_ref)) {
            stop("'max_cells_ref' must be a positive integer.")
        }
    }

    if (!is.null(max_cells_query)) {
        if (!is.numeric(max_cells_query) || max_cells_query <= 0 || max_cells_query != as.integer(max_cells_query)) {
            stop("'max_cells_query' must be a positive integer.")
        }
    }

    # Helper function to plot scores (original functionality with added downsampling)
    .plotScores <- function(x, cell_types,
                            sir_subset,
                            lower_facet,
                            diagonal_facet,
                            upper_facet,
                            max_cells_ref,
                            max_cells_query) {

        # Get SIR projections data
        sir_projections <- na.omit(x[["sir_projections"]])

        # Get available cell types if not specified by user
        if (is.null(cell_types)) {
            cell_types <- unique(sir_projections[["cell_type"]])
        } else {
            # Check if specified cell types exist in the data
            if (!all(cell_types %in% sir_projections[["cell_type"]])) {
                stop("One or more specified cell types not found in the data")
            }
            # Filter projections to include only specified cell types
            sir_projections <- sir_projections[sir_projections[["cell_type"]] %in%
                                                   cell_types, ]
        }

        # Separate reference and query data
        ref_data <- sir_projections[sir_projections[["dataset"]] == "Reference", ]
        query_data <- sir_projections[sir_projections[["dataset"]] == "Query", ]

        # Downsample reference data if max_cells_ref is specified
        if(!is.null(max_cells_ref)){
            # Input validation for max_cells_ref
            if (!is.numeric(max_cells_ref) || max_cells_ref <= 0 || max_cells_ref != as.integer(max_cells_ref)) {
                stop("'max_cells_ref' must be a positive integer.")
            }

            if(nrow(ref_data) > max_cells_ref){
                # Stratified sampling by cell type
                ref_data_list <- list()
                for(ct in cell_types){
                    ct_data <- ref_data[ref_data[["cell_type"]] == ct, ]
                    if(nrow(ct_data) > 0){
                        # Calculate proportional allocation
                        n_cells_ct <- min(nrow(ct_data),
                                          max(1, round(max_cells_ref * nrow(ct_data) / nrow(ref_data))))
                        if(nrow(ct_data) > n_cells_ct){
                            sampled_indices <- sample(nrow(ct_data), n_cells_ct)
                            ct_data <- ct_data[sampled_indices, ]
                        }
                        ref_data_list[[ct]] <- ct_data
                    }
                }
                ref_data <- do.call(rbind, ref_data_list)
                rownames(ref_data) <- NULL
            }
        }

        # Downsample query data if max_cells_query is specified
        if(!is.null(max_cells_query)){
            # Input validation for max_cells_query
            if (!is.numeric(max_cells_query) || max_cells_query <= 0 || max_cells_query != as.integer(max_cells_query)) {
                stop("'max_cells_query' must be a positive integer.")
            }

            if(nrow(query_data) > max_cells_query){
                # Stratified sampling by cell type
                query_data_list <- list()
                for(ct in cell_types){
                    ct_data <- query_data[query_data[["cell_type"]] == ct, ]
                    if(nrow(ct_data) > 0){
                        # Calculate proportional allocation
                        n_cells_ct <- min(nrow(ct_data),
                                          max(1, round(max_cells_query * nrow(ct_data) / nrow(query_data))))
                        if(nrow(ct_data) > n_cells_ct){
                            sampled_indices <- sample(nrow(ct_data), n_cells_ct)
                            ct_data <- ct_data[sampled_indices, ]
                        }
                        query_data_list[[ct]] <- ct_data
                    }
                }
                query_data <- do.call(rbind, query_data_list)
                rownames(query_data) <- NULL
            }
        }

        # Combine the potentially downsampled reference and query data
        sir_projections <- rbind(ref_data, query_data)

        # Create SIR column names with variance explained
        plot_names <- paste0(
            "SIR", sir_subset, " (",
            sprintf("%.1f%%", x[["percent_var"]][sir_subset]), ")")

        # Create a new data frame with selected SIR components
        sir_df <- data.frame(matrix(0, nrow = nrow(sir_projections),
                                    ncol = length(sir_subset)))
        colnames(sir_df) <- plot_names

        for (i in 1:length(sir_subset)) {
            sir_df[, i] <- sir_projections[, sir_subset[i]]
        }

        # Create a cell type dataset column for coloring
        cell_type_dataset <- paste(sir_projections[["dataset"]],
                                   sir_projections[["cell_type"]], sep = " ")

        # Define the order of cell type and dataset combinations
        order_combinations <- paste(
            rep(c("Reference", "Query"),
                length(cell_types)),
            rep(sort(cell_types), each = 2))

        cell_type_dataset <- factor(cell_type_dataset,
                                    levels = order_combinations)

        # Generate colors for cell types and datasets
        cell_type_colors <- generateColors(order_combinations, paired = TRUE)

        # Add the cell type dataset and colors columns to the SIR data frame
        sir_df[["cell_type_dataset"]] <- cell_type_dataset

        # Create a simple plot to extract the legend using GGally::grab_legend
        legend_plot <- ggplot2::ggplot(sir_df,
                                       ggplot2::aes(x = sir_df[, 1],
                                                    y = sir_df[, 2])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["cell_type_dataset"]])) +
            ggplot2::scale_color_manual(values = cell_type_colors,
                                        name = "Cell Type") +
            ggplot2::theme(
                legend.position = "right",
                legend.box = "vertical",
                legend.key = ggplot2::element_rect(fill = "white"),
                legend.background = ggplot2::element_blank()
            ) +
            ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))

        # Scatterplot facet function
        .scatterFunc <- function(data, mapping, ...) {

            ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::geom_point(alpha = 0.5, size = 1,
                                    ggplot2::aes(color = cell_type_dataset)) +
                ggplot2::scale_color_manual(values = cell_type_colors) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(
                        color = "black",
                        fill = NA,
                        linewidth = 0.5))
        }

        # Contour facet function
        .smoothContourFunc <- function(data, mapping, ...) {

            x_name <- rlang::as_name(mapping[["x"]])
            y_name <- rlang::as_name(mapping[["y"]])

            # Start with empty plot
            p <- ggplot2::ggplot() +
                ggplot2::theme_minimal()

            # Process each cell type separately
            for (ct in unique(data[["cell_type_dataset"]])) {
                subset_data <- data[data[["cell_type_dataset"]] == ct, ]

                # Skip if too few points
                if (nrow(subset_data) < 10) next

                # Use larger adjustment factor for smoother contours
                adjust_factor <- 1.5

                # Add smoother contours with fewer levels
                p <- p + ggplot2::stat_density_2d(
                    data = subset_data,
                    mapping = ggplot2::aes(
                        x = !!rlang::sym(x_name),
                        y = !!rlang::sym(y_name)
                    ),
                    contour = TRUE,
                    adjust = adjust_factor,
                    bins = 5,
                    color = cell_type_colors[which(
                        levels(data[["cell_type_dataset"]]) == ct)],
                    linewidth = 0.5,
                    na.rm = TRUE
                )
            }

            # Apply theme elements
            p + ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(
                        color = "black", fill = NA, linewidth = 0.5),
                    legend.position = "none",
                    axis.text = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_blank(),
                    axis.title = ggplot2::element_blank()
                )
        }

        # Ellipse facet function
        .robustEllipseFunc <- function(data, mapping, ...) {

            # Function to calculate robust ellipses through bootstrapping
            createEllipse <- function(d) {
                if (nrow(d) < 10) return(NULL)

                x_var <- rlang::as_name(mapping[["x"]])
                y_var <- rlang::as_name(mapping[["y"]])

                # Extract variables
                x <- d[[x_var]]
                y <- d[[y_var]]

                # Calculate robust center and covariance
                cov_mat <- cov(cbind(x, y), use = "pairwise.complete.obs")
                center <- c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))

                # Calculate ellipse points
                theta <- seq(0, 2 * pi, length.out = 100)
                ellipse <- MASS::cov.trob(cbind(x, y))

                # Get ellipse coordinates
                ev <- eigen(ellipse[["cov"]])
                a <- sqrt(ev[["values"]][1]) * 2.45
                b <- sqrt(ev[["values"]][2]) * 2.45

                # Create ellipse coordinates
                angle <- atan2(ev[["vectors"]][2,1], ev[["vectors"]][1,1])
                ellipse_x <- center[1] + a * cos(theta) * cos(angle) -
                    b * sin(theta) * sin(angle)
                ellipse_y <- center[2] + a * cos(theta) * sin(angle) +
                    b * sin(theta) * cos(angle)

                return(data.frame(x = ellipse_x, y = ellipse_y))
            }

            p <- ggplot2::ggplot(data, mapping)

            # Split by cell type and create robust ellipses
            for (ct in unique(data[["cell_type_dataset"]])) {
                subset_data <- data[data[["cell_type_dataset"]] == ct,]
                ellipse_data <- createEllipse(subset_data)

                if (!is.null(ellipse_data)) {
                    p <- p + ggplot2::geom_path(
                        data = ellipse_data,
                        ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                        color = cell_type_colors[
                            which(levels(data[["cell_type_dataset"]]) == ct)],
                        linewidth = 0.7
                    )
                }
            }

            p + ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(
                        color = "black", fill = NA,
                        linewidth = 0.5),
                    legend.position = "none",
                    axis.text = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_blank(),
                    axis.title = ggplot2::element_blank()
                )
        }

        # Blank facet function
        .blankFunc <- function(data, mapping, ...) {

            ggplot2::ggplot() +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(
                        color = "black", fill = NA,
                        linewidth = 0.5),
                    legend.position = "none",
                    axis.text = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_blank(),
                    axis.title = ggplot2::element_blank(),
                    panel.grid = ggplot2::element_blank()
                )
        }

        # Ridge diagonal facet
        .ridgeFunc <- function(data, mapping, ...) {

            # Get current mapping info
            x_var_name <- rlang::as_name(mapping[["x"]])

            # Create a long-format data frame for ridge plot
            plot_data <- data.frame(
                value = data[[x_var_name]],
                group = data[["cell_type_dataset"]]
            )

            # Create ridge plot with ggridges
            suppressMessages({
                p <- ggplot2::ggplot(plot_data,
                                     ggplot2::aes(x = .data[["value"]],
                                                  y = .data[["group"]],
                                                  fill = .data[["group"]])) +
                    ggridges::geom_density_ridges(
                        alpha = 0.7,
                        scale = 2,
                        rel_min_height = 0.01,
                        quantile_lines = FALSE
                    ) +
                    ggplot2::scale_fill_manual(values = cell_type_colors) +
                    ggplot2::scale_y_discrete(
                        limits = rev(levels(cell_type_dataset))) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(
                        panel.border = ggplot2::element_rect(color = "black",
                                                             fill = NA,
                                                             linewidth = 0.5),
                        axis.title = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        legend.position = "none"
                    )
            })

            return(p)
        }

        # Density diagonal facet
        .densityFunc <- function(data, mapping, ...) {

            # Get current mapping info
            x_var_name <- rlang::as_name(mapping[["x"]])

            ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::geom_density(ggplot2::aes(fill = cell_type_dataset,
                                                   color = cell_type_dataset),
                                      alpha = 0.5) +
                ggplot2::scale_fill_manual(values = cell_type_colors) +
                ggplot2::scale_color_manual(values = cell_type_colors) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(color = "black",
                                                         fill = NA,
                                                         linewidth = 0.5),
                    axis.title.y = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank()
                )
        }

        # Boxplot diagonal facet
        .boxplotFunc <- function(data, mapping, ...) {

            # Extract the x variable name for the boxplot title
            x_name <- rlang::as_name(mapping[["x"]])

            # Create a long-format data frame for horizontal boxplot
            plot_data <- data.frame(
                value = data[[x_name]],
                group = data[["cell_type_dataset"]]
            )

            # Create horizontal boxplot
            ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["value"]],
                                                    y = .data[["group"]],
                                                    fill = .data[["group"]])) +
                ggplot2::geom_boxplot(alpha = 0.7,
                                      outlier.size = 0.5,
                                      width = 0.6) +
                ggplot2::scale_fill_manual(values = cell_type_colors) +
                ggplot2::scale_y_discrete(limits = rev(
                    levels(cell_type_dataset))) +
                ggplot2::labs(x = "", y = "") +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(color = "black",
                                                         fill = NA,
                                                         linewidth = 0.5),
                    axis.title.x = ggplot2::element_blank(),
                    legend.position = "none",
                    axis.title.y = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank()
                )
        }

        # Select the lower facet plot type based on user input
        if (lower_facet == "scatter") {
            lower_plot <- .scatterFunc
        } else if (lower_facet == "contour") {
            lower_plot <- .smoothContourFunc
        } else if (lower_facet == "ellipse") {
            lower_plot <- .robustEllipseFunc
        } else if (lower_facet == "blank") {
            lower_plot <- .blankFunc
        }

        # Select the diagonal plot type based on user input
        if (diagonal_facet == "density") {
            diag_plot <- .densityFunc
        } else if (diagonal_facet == "boxplot") {
            diag_plot <- .boxplotFunc
        } else if (diagonal_facet == "ridge") {
            diag_plot <- .ridgeFunc
        }

        # Select the upper facet plot type based on user input
        if (upper_facet == "blank") {
            upper_plot <- .blankFunc
        } else if (upper_facet == "scatter") {
            upper_plot <- .scatterFunc
        } else if (upper_facet == "contour") {
            upper_plot <- .smoothContourFunc
        } else if (upper_facet == "ellipse") {
            upper_plot <- .robustEllipseFunc
        }

        # Create pairs plot using GGally
        plot_obj <- suppressMessages(
            GGally::ggpairs(
                sir_df,
                columns = seq_len(length(sir_subset)),
                mapping = ggplot2::aes(color = cell_type_dataset),
                lower = list(continuous = lower_plot),
                upper = list(continuous = upper_plot),
                diag = list(continuous = diag_plot),
                progress = FALSE,
                legend = GGally::grab_legend(legend_plot)
            )
        )

        # Add black frame around facet titles
        plot_obj <- plot_obj +
            ggplot2::theme(
                strip.background = ggplot2::element_rect(
                    fill = "white", color = "black", linewidth = 0.5),
                strip.text = ggplot2::element_text(color = "black")
            )

        # Return the plot
        return(plot_obj)
    }

    # Helper function to plot loadings (unchanged)
    .plotLoadings <- function(x, sir_subset, n_top) {

        # Get rotation matrix (loadings)
        rotation_mat <- x[["rotation_mat"]]

        # Check if rotation matrix exists
        if (is.null(rotation_mat)) {
            stop("No rotation matrix found in the SIR object.")
        }

        # Get variable names
        var_names <- rownames(rotation_mat)
        if (is.null(var_names)) {
            var_names <- paste0("Var_", seq_len(nrow(rotation_mat)))
        }

        # Create a combined data frame for all SIR components
        all_loadings_list <- list()

        for (i in seq_along(sir_subset)) {
            sir_comp <- sir_subset[i]

            # Get loadings for this component
            loadings <- rotation_mat[, sir_comp]

            # Create data frame with variable names and loadings
            loading_df <- data.frame(
                variable = var_names,
                loading = loadings,
                abs_loading = abs(loadings),
                sir_component = paste0("SIR", sir_comp, " (",
                                       sprintf("%.1f%%",
                                               x[["percent_var"]][sir_comp]),
                                       ")"),
                sir_order = i,
                stringsAsFactors = FALSE
            )

            # Sort by absolute loading and take n_top variables
            loading_df <- loading_df[order(loading_df[["abs_loading"]],
                                           decreasing = TRUE), ]
            loading_df <- loading_df[seq_len(min(n_top,
                                                 nrow(loading_df))), ]

            # Order variables for plotting (highest absolute loading at n_top of each facet)
            loading_df <- loading_df[order(loading_df[["abs_loading"]]), ]
            loading_df[["variable_facet"]] <- factor(loading_df[["variable"]],
                                                     levels = loading_df[["variable"]])

            # Create color based on sign of loading
            loading_df[["color"]] <- ifelse(loading_df[["loading"]] >= 0,
                                            "Positive", "Negative")

            all_loadings_list[[i]] <- loading_df
        }

        # Combine all data frames
        all_loadings <- do.call(rbind, all_loadings_list)

        # Create factor for SIR components to control facet order
        sir_component_levels <- unique(all_loadings[["sir_component"]][
            order(all_loadings[["sir_order"]])])
        all_loadings[["sir_component"]] <- factor(
            all_loadings[["sir_component"]],
            levels = sir_component_levels)

        # Create the faceted plot
        p <- ggplot2::ggplot(all_loadings, ggplot2::aes(
            x = .data[["loading"]],
            y = .data[["variable_facet"]],
            fill = .data[["color"]])) +
            ggplot2::geom_col(alpha = 0.8, width = 0.6) +
            ggplot2::facet_wrap(~ sir_component, scales = "free_y",
                                ncol = length(sir_subset)) +
            ggplot2::scale_fill_manual(values = c("Positive" = "#2166ac",
                                                  "Negative" = "#d6604d"),
                                       name = "Loading") +
            ggplot2::labs(
                title = "SIR Component Loadings",
                x = "Loading Value",
                y = "Variables"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(color = "black",
                                                     fill = NA,
                                                     linewidth = 0.5),
                strip.background = ggplot2::element_rect(fill = "white",
                                                         color = "black",
                                                         linewidth = 0.5),
                strip.text = ggplot2::element_text(size = 10,
                                                   face = "bold"),
                plot.title = ggplot2::element_text(hjust = 0.5,
                                                   size = 14,
                                                   face = "bold"),
                axis.text.y = ggplot2::element_text(size = 8),
                legend.position = "bottom"
            ) +
            ggplot2::geom_vline(xintercept = 0,
                                linetype = "dashed",
                                color = "gray50",
                                alpha = 0.7)

        return(p)
    }


    # Branch based on plot type
    if (plot_type == "scores") {
        return(.plotScores(x,
                           cell_types,
                           sir_subset,
                           lower_facet,
                           diagonal_facet,
                           upper_facet,
                           max_cells_ref,
                           max_cells_query))
    } else {
        return(.plotLoadings(x,
                             sir_subset,
                             n_top))
    }
}
