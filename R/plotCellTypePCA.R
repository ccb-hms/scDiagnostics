#' @title Plot Principal Components for Different Cell Types
#'
#' @description
#' This function plots the principal components for different cell types in the query and reference datasets.
#'
#' @details
#' This function projects the query dataset onto the principal component space of the reference dataset and then plots the
#' specified principal components for the specified cell types.
#' It uses the `projectPCA` function to perform the projection and \code{GGally} to create the pairs plot.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is 1:5.
#' @param lower_facet Type of plot to use for the lower panels. Either "scatter" (default), "contour", "ellipse", or "blank".
#' @param diagonal_facet Type of plot to use for the diagonal panels. Either "ridge" (default), "density", or "boxplot".
#' @param upper_facet Type of plot to use for the upper panels. Either "blank" (default), "scatter", "contour", or "ellipse".
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A ggmatrix object representing a pairs plot of specified principal components for the given cell types and datasets.
#'
#' @export
#'
# Function to plot PC for different cell types
plotCellTypePCA <- function(query_data,
                            reference_data,
                            query_cell_type_col,
                            ref_cell_type_col,
                            cell_types = NULL,
                            pc_subset = 1:5,
                            assay_name = "logcounts",
                            lower_facet = c("scatter", "contour", "ellipse", "blank"),
                            diagonal_facet = c("ridge", "density", "boxplot"),
                            upper_facet = c("blank", "scatter", "contour", "ellipse"),
                            max_cells = 2500){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)

    # Match diagonal_facet and upper_facet arguments
    lower_facet <- match.arg(lower_facet)
    diagonal_facet <- match.arg(diagonal_facet)
    upper_facet <- match.arg(upper_facet)

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             pc_subset = pc_subset,
                             assay_name = assay_name)
    pca_output <- pca_output[pca_output[["cell_type"]] %in% cell_types,]

    # Create PC column names with variance explained
    plot_names <- paste0(
        "PC", pc_subset, " (",
        sprintf("%.1f%%", attributes(
            reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset]), ")")

    # Create a new data frame with selected PCs
    pc_df <- data.frame(matrix(0, nrow = nrow(pca_output),
                               ncol = length(pc_subset)))
    colnames(pc_df) <- plot_names

    for (i in 1:length(pc_subset)) {
        pc_df[, i] <- pca_output[, paste0("PC", pc_subset[i])]
    }

    # Create a cell type dataset column for coloring
    cell_type_dataset <- paste(pca_output[["dataset"]],
                               pca_output[["cell_type"]], sep = " ")

    # Define the order of cell type and dataset combinations
    order_combinations <- paste(
        rep(c("Reference", "Query"),
            length(cell_types)),
        rep(sort(cell_types), each = 2))

    cell_type_dataset <- factor(cell_type_dataset,
                                levels = order_combinations)

    # Generate colors for cell types and datasets
    cell_type_colors <- generateColors(order_combinations, paired = TRUE)

    # Add the cell type dataset column to the PC data frame
    pc_df[["cell_type_dataset"]] <- cell_type_dataset

    # Create a simple plot to extract the legend using GGally::grab_legend
    legend_plot <- ggplot2::ggplot(pc_df,
                                   ggplot2::aes(x = pc_df[,1],
                                                y = pc_df[,2])) +
        ggplot2::geom_point(ggplot2::aes(color = cell_type_dataset)) +
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

    # Select the upper facet plot type based on user input
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
            pc_df,
            columns = seq_len(length(pc_subset)),
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
