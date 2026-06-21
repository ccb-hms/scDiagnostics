#' @title Create Visualization Plots From \code{detectAnomaly} Object
#'
#' @description
#' This S3 plot method generates visualizations for anomaly detection results.
#' If PCA was used (\code{pc_subset} is numeric), it generates faceted scatter plots for the principal components.
#' If highly variable genes (HVGs) were used (\code{pc_subset} is NULL), it generates a ComplexHeatmap of the HVGs.
#'
#' @details
#' **PCA Scatter Plots:**
#' Extracts the specified PCs and generates a `GGally` pairs plot. Lower facets show scatter plots
#' with a background gradient representing anomaly scores. Diagonal facets show distributions, and
#' upper facets show contours/ellipses separated by anomaly status.
#'
#' **HVG Heatmaps:**
#' Extracts the highly variable genes (HVGs) and generates a `ComplexHeatmap`. Cells are ordered by
#' Dataset (Query vs Reference) and by Anomaly Status (Anomalous vs Non-Anomalous). Gene expression
#' is Z-score scaled across cells for optimal visual contrast.
#'
#' @param x A list object containing the anomaly detection results from the \code{detectAnomaly} function.
#' @param cell_type A character string specifying the cell type for which the plots should be generated.
#' If NULL, the "Combined" cell type will be plotted. Default is NULL.
#' @param pc_subset A numeric vector specifying the indices of the PCs to be included in the plots.
#' If NULL, all PCs in \code{reference_mat_subset} will be included. Ignored if HVGs were used.
#' @param data_type A character string specifying whether to plot the "query" data, "reference" data,
#' or "both". Note: "both" is only supported for HVG Heatmaps. Default is "query".
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 500.
#' @param upper_facet Either "blank" (default), "contour", or "ellipse" for the upper facet plots (PCA only).
#' @param diagonal_facet Either "density" (default), "ridge", "boxplot" or "blank" for the diagonal plots (PCA only).
#' @param max_cells_ref Maximum number of reference cells to include in the plot. If NULL, all are plotted. Default is NULL.
#' @param max_cells_query Maximum number of query cells to include in the plot. If NULL, all are plotted. Default is NULL.
#' @param draw_plot Logical indicating whether to draw the plot immediately (TRUE) or return
#' the undrawn plot object (FALSE). For heatmaps, FALSE returns a ComplexHeatmap object. Default is TRUE.
#' @param ... Additional arguments passed to the `isolation.forest` function (PCA) or `ComplexHeatmap::Heatmap` function.
#'
#' @return Returns a \code{GGally::ggpairs} object for PCA data, or a \code{ComplexHeatmap} object for HVG data.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{detectAnomaly}}
#'
#' @rdname detectAnomaly
#'
# Function to plot anomaly regions or heatmaps
plot.detectAnomalyObject <- function(x,
                                     cell_type = NULL,
                                     pc_subset = NULL,
                                     data_type = c("query", "reference", "both"),
                                     n_tree = 500,
                                     upper_facet = c("blank", "contour", "ellipse"),
                                     diagonal_facet = c("density", "ridge", "boxplot", "blank"),
                                     max_cells_ref = NULL,
                                     max_cells_query = NULL,
                                     draw_plot = TRUE,
                                     ...) {

    # Check input for cell type
    if(is.null(cell_type)){
        cell_type <- "Combined"
    } else {
        if(!(cell_type %in% names(x)))
            stop("'cell_type' is not available in 'x'.")
    }

    # Determine if this is PCA data or HVG data
    is_pca <- "var_explained" %in% names(x[[cell_type]])

    # Match arguments
    data_type <- match.arg(data_type)
    upper_facet <- match.arg(upper_facet)
    diagonal_facet <- match.arg(diagonal_facet)

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

    # __________________________________
    # PATH A: HVG COMPLEX HEATMAP LOGIC
    # __________________________________
    if (!is_pca) {

        if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
            stop("Packages 'ComplexHeatmap' and 'circlize' are required to plot HVG heatmaps. Please install them.")
        }

        # Check query availability
        if(is.null(x[[cell_type]][["query_mat_subset"]]) && data_type %in% c("query", "both")){
            stop("There is no query data available in the 'detectAnomaly' object.")
        }

        # Extract matrices and anomaly statuses
        ref_mat <- x[[cell_type]][["reference_mat_subset"]]
        ref_anomaly <- x[[cell_type]][["reference_anomaly"]]

        query_mat <- x[[cell_type]][["query_mat_subset"]]
        query_anomaly <- x[[cell_type]][["query_anomaly"]]

        # Downsample if requested
        if(!is.null(max_cells_ref) && nrow(ref_mat) > max_cells_ref){
            idx <- sample(nrow(ref_mat), max_cells_ref)
            ref_mat <- ref_mat[idx, , drop = FALSE]
            ref_anomaly <- ref_anomaly[idx]
        }
        if(!is.null(max_cells_query) && !is.null(query_mat) && nrow(query_mat) > max_cells_query){
            idx <- sample(nrow(query_mat), max_cells_query)
            query_mat <- query_mat[idx, , drop = FALSE]
            query_anomaly <- query_anomaly[idx]
        }

        # Combine data based on data_type
        if (data_type == "both") {
            mat <- rbind(ref_mat, query_mat)
            dataset <- c(rep("Reference", nrow(ref_mat)), rep("Query", nrow(query_mat)))

            # Force reference anomalies to FALSE so we don't distract the user
            anomaly <- c(rep(FALSE, length(ref_anomaly)), query_anomaly)

            cell_names <- c(paste0("Ref_", seq_len(nrow(ref_mat))), paste0("Query_", seq_len(nrow(query_mat))))
        } else if (data_type == "query") {
            mat <- query_mat
            dataset <- rep("Query", nrow(query_mat))
            anomaly <- query_anomaly
            cell_names <- paste0("Query_", seq_len(nrow(query_mat)))
        } else {
            mat <- ref_mat
            dataset <- rep("Reference", nrow(ref_mat))
            anomaly <- ref_anomaly
            cell_names <- paste0("Ref_", seq_len(nrow(ref_mat)))
        }

        rownames(mat) <- cell_names
        plot_mat <- t(mat) # Rows = Genes, Columns = Cells

        # Z-score scaling by gene (row)
        plot_mat <- t(scale(t(plot_mat)))
        plot_mat[is.na(plot_mat)] <- 0 # Handle zero-variance genes
        plot_mat[plot_mat > 2] <- 2    # Cap extremes for better color contrast
        plot_mat[plot_mat < -2] <- -2

        # Determine anomaly labels
        anomaly_labels <- ifelse(anomaly, "Anomalous", "Non-Anomalous")

        # Create annotation dataframe
        annotation_col <- data.frame(
            Dataset = factor(dataset, levels = c("Query", "Reference")),
            Status = factor(anomaly_labels, levels = c("Anomalous", "Non-Anomalous"))
        )
        rownames(annotation_col) <- cell_names

        # Order columns clearly (Query Anomalous -> Query Normal -> Reference Anomalous -> Reference Normal)
        order_idx <- order(annotation_col$Dataset, annotation_col$Status)
        plot_mat <- plot_mat[, order_idx, drop = FALSE]
        annotation_col <- annotation_col[order_idx, , drop = FALSE]

        # Define colors (Matching calculateGeneShiftsObject styling)
        dataset_colors <- c("Query" = "#B565D8", "Reference" = "#5A9BD8")
        anomaly_colors <- c("Non-Anomalous" = "#9E9E9E", "Anomalous" = "#D2314C")

        zscore_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

        # Create Top Annotation
        top_ha <- ComplexHeatmap::HeatmapAnnotation(
            Status = annotation_col$Status,
            Dataset = annotation_col$Dataset,
            col = list(Status = anomaly_colors, Dataset = dataset_colors),
            show_legend = TRUE,
            show_annotation_name = TRUE,
            annotation_name_side = "left"
        )

        # Build Heatmap
        ht <- ComplexHeatmap::Heatmap(
            matrix = plot_mat,
            name = "Z-Score",
            col = zscore_col_fun,
            heatmap_legend_param = list(
                title = "Z-Score",
                at = c(-2, -1, 0, 1, 2),
                labels = c("-2", "-1", "0", "1", "2")
            ),
            # Only show gene names if there aren't too many of them
            show_row_names = nrow(plot_mat) <= 100,
            row_names_gp = if (requireNamespace("grid", quietly = TRUE)) grid::gpar(fontsize = 8) else NULL,
            cluster_rows = TRUE,
            show_column_names = FALSE,
            cluster_columns = FALSE, # Columns are strictly ordered by Dataset & Status
            column_order = colnames(plot_mat),
            top_annotation = top_ha,
            border = TRUE,
            column_title = paste0("HVG Expression Heatmap: ", cell_type),
            use_raster = ncol(plot_mat) > 1000 || nrow(plot_mat) > 1000,
            raster_quality = 2,
            ...
        )

        if (draw_plot) {
            return(ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE))
        } else {
            return(ht)
        }
    }

    # _______________________________
    # PATH B: PCA SCATTER PLOT LOGIC
    # _______________________________

    if (data_type == "both") {
        stop("data_type = 'both' is only supported for HVG heatmaps. For PCA plots, please specify 'query' or 'reference'.")
    }

    # Check input for pc_subset
    if(!is.null(pc_subset)){
        if(!all(pc_subset %in% seq_len(length(x[[cell_type]][["var_explained"]]))))
            stop("'pc_subset' is out of range.")
    } else {
        pc_subset <- seq_len(length(x[[cell_type]][["var_explained"]]))
    }

    # Filter and prepare data based on data type
    if(is.null(x[[cell_type]][["query_mat_subset"]]) && data_type == "query"){
        stop("There is no query data available in the 'detectAnomaly' object.")
    } else {
        if(data_type == "query"){
            data_subset <- x[[cell_type]][["query_mat_subset"]][, pc_subset, drop = FALSE]
            anomaly <- x[[cell_type]][["query_anomaly"]]

            if(!is.null(max_cells_query) && nrow(data_subset) > max_cells_query){
                sampled_indices <- sample(nrow(data_subset), max_cells_query)
                data_subset <- data_subset[sampled_indices, , drop = FALSE]
                anomaly <- anomaly[sampled_indices]
            }
        } else if(data_type == "reference"){
            data_subset <- x[[cell_type]][["reference_mat_subset"]][, pc_subset, drop = FALSE]
            anomaly <- x[[cell_type]][["reference_anomaly"]]

            if(!is.null(max_cells_ref) && nrow(data_subset) > max_cells_ref){
                sampled_indices <- sample(nrow(data_subset), max_cells_ref)
                data_subset <- data_subset[sampled_indices, , drop = FALSE]
                anomaly <- anomaly[sampled_indices]
            }
        }
    }

    # Add variance explained to column names
    colnames(data_subset) <- paste0(
        "PC", pc_subset,
        " (", sprintf("%.1f%%", x[[cell_type]][["var_explained"]][pc_subset]), ")")

    # Create a data frame with PC values and anomaly info
    pc_df <- data.frame(data_subset)
    colnames(pc_df) <- colnames(data_subset)
    pc_df[["anomaly"]] <- factor(anomaly)

    # Colors for anomalous and non-anomalous data
    anomaly_colors <- c("TRUE" = "red", "FALSE" = "black")
    anomaly_fill_colors <- c("TRUE" = "#FFB6B6", "FALSE" = "#DDDDDD")

    # Train isolation forest for each PC combination
    isolation_forests <- list()
    for(i in seq_along(pc_subset)){
        for(j in seq_along(pc_subset)){
            if(i < j){
                pc_i <- paste0("PC", pc_subset[i])
                pc_j <- paste0("PC", pc_subset[j])

                train_data <- x[[cell_type]][["reference_mat_subset"]][, c(pc_i, pc_j)]

                isolation_forests[[paste(pc_i, pc_j, sep="-")]] <-
                    isotree::isolation.forest(train_data, ntree = n_tree, ...)
            }
        }
    }

    # --- Plotting helper functions ---
    .anomalyScatterFunc <- function(data, mapping, ...){
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])

        pc_x <- as.numeric(sub("PC([0-9]+).*", "\\1", x_name))
        pc_y <- as.numeric(sub("PC([0-9]+).*", "\\1", y_name))

        if_key <- paste0("PC", min(pc_x, pc_y), "-PC", max(pc_x, pc_y))
        if_model <- isolation_forests[[if_key]]

        x_range <- range(data[[x_name]])
        y_range <- range(data[[y_name]])
        x_buffer <- 0.1 * diff(x_range)
        y_buffer <- 0.1 * diff(y_range)

        x_seq <- seq(x_range[1] - x_buffer, x_range[2] + x_buffer, length.out = 150)
        y_seq <- seq(y_range[1] - y_buffer, y_range[2] + y_buffer, length.out = 150)
        grid <- expand.grid(x = x_seq, y = y_seq)

        x_step <- diff(x_seq)[1]
        y_step <- diff(y_seq)[1]

        grid_data <- data.frame(grid)
        colnames(grid_data) <- c(paste0("PC", pc_x), paste0("PC", pc_y))

        scores <- predict(if_model, newdata = grid_data, type = "score")

        background_data <- data.frame(
            x_value = grid[["x"]],
            y_value = grid[["y"]],
            anomaly_score = scores,
            width = x_step,
            height = y_step
        )

        gradient_colors <- c("#2CA25F", "#FFFFE5", "#FF8585")
        stops <- c(0, 0.55, 1)

        p <- ggplot2::ggplot() +
            ggplot2::geom_tile(
                data = background_data,
                ggplot2::aes(x = .data[["x_value"]], y = .data[["y_value"]],
                             fill = .data[["anomaly_score"]], width = .data[["width"]], height = .data[["height"]])
            ) +
            ggplot2::scale_fill_gradientn(
                colors = gradient_colors, values = scales::rescale(stops), limits = c(0, 1), guide = "none"
            ) +
            ggplot2::geom_point(
                data = data,
                ggplot2::aes(x = !!rlang::sym(x_name), y = !!rlang::sym(y_name), color = anomaly),
                shape = 16, size = 1.5, alpha = 0.45
            ) +
            ggplot2::scale_color_manual(values = anomaly_colors, name = "Anomalous") +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = c(0, 0)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(),
                legend.position = "none"
            )
        return(p)
    }

    .densityDiagFunc <- function(data, mapping, ...){
        x_name <- rlang::as_name(mapping[["x"]])
        ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_name))) +
            ggplot2::geom_density(ggplot2::aes(fill = anomaly, color = anomaly), alpha = 0.3, linewidth = 0.8) +
            ggplot2::scale_color_manual(values = anomaly_colors) +
            ggplot2::scale_fill_manual(values = anomaly_fill_colors) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            )
    }

    .ridgeDiagFunc <- function(data, mapping, ...){
        x_name <- rlang::as_name(mapping[["x"]])
        plot_data <- data.frame(value = data[[x_name]], anomaly = data[["anomaly"]])
        suppressMessages({
            p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["value"]], y = .data[["anomaly"]],
                                                         fill = .data[["anomaly"]], color = .data[["anomaly"]])) +
                ggridges::geom_density_ridges(alpha = 0.6, scale = 1.8, rel_min_height = 0.01, quantile_lines = FALSE) +
                ggplot2::scale_fill_manual(values = anomaly_fill_colors) +
                ggplot2::scale_color_manual(values = anomaly_colors) +
                ggplot2::scale_y_discrete(limits = rev(levels(data[["anomaly"]]))) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                    axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(), legend.position = "none",
                    panel.grid = ggplot2::element_blank()
                )
        })
        return(p)
    }

    .boxplotDiagFunc <- function(data, mapping, ...){
        x_name <- rlang::as_name(mapping[["x"]])
        plot_data <- data.frame(value = data[[x_name]], anomaly = data[["anomaly"]])
        ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["value"]], y = .data[["anomaly"]], fill = .data[["anomaly"]])) +
            ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6, color = "black") +
            ggplot2::scale_fill_manual(values = anomaly_fill_colors) +
            ggplot2::scale_y_discrete(limits = rev(levels(data[["anomaly"]]))) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            )
    }

    .blankFunc <- function(data, mapping, ...){
        ggplot2::ggplot() + ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()
            )
    }

    .contourFunc <- function(data, mapping, ...){
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])
        p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_name), y = !!rlang::sym(y_name))) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none"
            )
        data_normal <- data[data[["anomaly"]] == "FALSE",]
        if(nrow(data_normal) >= 10) {
            p <- p + ggplot2::stat_density_2d(data = data_normal, ggplot2::aes(x = !!rlang::sym(x_name), y = !!rlang::sym(y_name)),
                                              color = anomaly_colors[["FALSE"]], linewidth = 0.5, contour = TRUE, bins = 5)
        }
        data_anomaly <- data[data[["anomaly"]] == "TRUE",]
        if(nrow(data_anomaly) >= 10) {
            p <- p + ggplot2::stat_density_2d(data = data_anomaly, ggplot2::aes(x = !!rlang::sym(x_name), y = !!rlang::sym(y_name)),
                                              color = anomaly_colors[["TRUE"]], linewidth = 0.5, contour = TRUE, bins = 5)
        }
        return(p)
    }

    .ellipseFunc <- function(data, mapping, ...){
        x_name <- rlang::as_name(mapping[["x"]])
        y_name <- rlang::as_name(mapping[["y"]])
        createEllipse <- function(d) {
            if (nrow(d) < 10) return(NULL)
            x <- d[[x_name]]; y <- d[[y_name]]
            cov_mat <- cov(cbind(x, y), use = "pairwise.complete.obs")
            center <- c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
            ellipse <- MASS::cov.trob(cbind(x, y))
            ev <- eigen(ellipse[["cov"]])
            a <- sqrt(ev[["values"]][1]) * 2.45
            b <- sqrt(ev[["values"]][2]) * 2.45
            # BUG FIX: Use 'vectors' instead of 'values' for correct dimension
            angle <- atan2(ev[["vectors"]][2,1], ev[["vectors"]][1,1])
            theta <- seq(0, 2 * pi, length.out = 100)
            ellipse_x <- center[1] + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
            ellipse_y <- center[2] + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)
            return(data.frame(x = ellipse_x, y = ellipse_y))
        }
        p <- ggplot2::ggplot() + ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none"
            )
        data_normal <- data[data[["anomaly"]] == "FALSE",]
        ellipse_normal <- createEllipse(data_normal)
        if (!is.null(ellipse_normal)) {
            p <- p + ggplot2::geom_path(data = ellipse_normal, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                                        color = anomaly_colors[["FALSE"]], linewidth = 0.7)
        }
        data_anomaly <- data[data[["anomaly"]] == "TRUE",]
        ellipse_anomaly <- createEllipse(data_anomaly)
        if (!is.null(ellipse_anomaly)) {
            p <- p + ggplot2::geom_path(data = ellipse_anomaly, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                                        color = anomaly_colors[["TRUE"]], linewidth = 0.7)
        }
        return(p)
    }

    if (diagonal_facet == "density") { diag_func <- .densityDiagFunc }
    else if (diagonal_facet == "ridge") { diag_func <- .ridgeDiagFunc }
    else if (diagonal_facet == "boxplot") { diag_func <- .boxplotDiagFunc }
    else if (diagonal_facet == "blank") { diag_func <- .blankFunc }

    if (upper_facet == "blank") { upper_func <- .blankFunc }
    else if (upper_facet == "contour") { upper_func <- .contourFunc }
    else if (upper_facet == "ellipse") { upper_func <- .ellipseFunc }

    legend_plot <- ggplot2::ggplot(pc_df, ggplot2::aes(x = pc_df[,1], y = pc_df[,2])) +
        ggplot2::geom_point(ggplot2::aes(color = anomaly)) +
        ggplot2::scale_color_manual(values = anomaly_colors, name = "Anomalous") +
        ggplot2::theme(legend.position = "right", legend.box = "vertical", legend.key = ggplot2::element_rect(fill = "white"))

    plot_obj <- suppressMessages(
        GGally::ggpairs(
            pc_df, columns = seq_len(length(pc_subset)), mapping = ggplot2::aes(color = anomaly),
            lower = list(continuous = .anomalyScatterFunc), diag = list(continuous = diag_func),
            upper = list(continuous = upper_func), progress = FALSE, legend = GGally::grab_legend(legend_plot)
        )
    )

    plot_obj <- plot_obj +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5),
            strip.text = ggplot2::element_text(color = "black"), plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::ggtitle(paste0("Isolation Forest Anomaly Plot: ", cell_type))

    if (draw_plot) {
        print(plot_obj)
        return(invisible(plot_obj))
    } else {
        return(plot_obj)
    }
}
