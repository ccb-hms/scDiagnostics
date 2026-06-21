#' @title Plot PCA Reconstruction Errors
#'
#' @description
#' This S3 plot method generates visualizations for PCA reconstruction errors (out-of-distribution anomalies).
#' It displays the distribution of errors for the query and/or reference datasets, highlights the
#' calculated anomaly threshold, and explicitly identifies anomalous cells. It can also generate
#' heatmaps of the highly variable genes driving the local PCA space.
#'
#' @details
#' The function extracts the reconstruction errors from the given object and generates a visualization.
#'
#' Four \code{plot_type} options are available:
#' \itemize{
#'   \item \code{"violin"} (Default): Shows a violin plot of the error distribution overlaid with individual cell
#'         points (jittered) colored by anomaly status. The violin is trimmed to the data range to maintain
#'         statistical validity (preventing density estimation below zero).
#'   \item \code{"boxplot"}: Shows a standard boxplot overlaid with jittered points.
#'   \item \code{"ridge"}: Shows ridgeline plots separating the datasets vertically. A vertical red dashed line
#'         marks the threshold. Best for visualizing the density of the non-anomalous cells versus the long tail of anomalies.
#'   \item \code{"heatmap"}: Generates a \code{ComplexHeatmap} showing the Z-score scaled expression of the
#'         highly variable genes used to construct the local PCA space. Cells are grouped by Dataset and Anomaly Status.
#' }
#'
#' @param x A list object of class \code{calculateReconstructionErrorObject} containing the results
#' from the \code{calculateReconstructionError} function.
#' @param cell_type A character string specifying the cell type for which the plots should be generated.
#' If NULL, defaults to "Combined" if available, otherwise plots the first available cell type. Default is NULL.
#' @param data_type A character string specifying whether to plot the "query" data, "reference" data,
#' or "both". Default is "both".
#' @param plot_type A character string specifying the type of visualization. Options are
#' \code{"violin"}, \code{"boxplot"}, \code{"ridge"}, or \code{"heatmap"}. Default is \code{"violin"}.
#' @param draw_plot Logical indicating whether to draw the plot immediately (TRUE) or return
#' the undrawn plot object (FALSE). For heatmaps, FALSE returns a ComplexHeatmap object. Default is FALSE.
#' @param ... Additional arguments passed to \code{ComplexHeatmap::Heatmap} when \code{plot_type = "heatmap"}.
#'
#' @return A \code{ggplot2} object for distribution plots, or a \code{ComplexHeatmap} object for heatmaps.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateReconstructionError}}
#'
#' @rdname calculateReconstructionError
#'
# Function to plot reconstruction errors
plot.calculateReconstructionErrorObject <- function(x,
                                                    cell_type = NULL,
                                                    data_type = c("both", "query", "reference"),
                                                    plot_type = c("violin", "boxplot", "ridge", "heatmap"),
                                                    draw_plot = FALSE,
                                                    ...) {

    # Match arguments
    data_type <- match.arg(data_type)
    plot_type <- match.arg(plot_type)

    # Smart default for cell_type
    if(is.null(cell_type)){
        if ("Combined" %in% names(x)) {
            cell_type <- "Combined"
        } else {
            cell_type <- names(x)[1] # Fallback to the first available cell type
            message("Defaulting to cell_type: '", cell_type, "'")
        }
    } else {
        if(!(cell_type %in% names(x))) {
            stop(paste0("Cell type '", cell_type, "' is not available in the provided object."))
        }
    }

    # ______________________
    # PATH A: HEATMAP LOGIC
    # ______________________

    if (plot_type == "heatmap") {

        if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
            stop("Packages 'ComplexHeatmap' and 'circlize' are required to plot heatmaps. Please install them.")
        }

        if(is.null(x[[cell_type]][["query_mat_subset"]]) && data_type %in% c("query", "both")){
            stop("There is no query data available in the object to plot.")
        }

        ref_mat <- x[[cell_type]][["reference_mat_subset"]]
        ref_anomaly <- x[[cell_type]][["reference_anomaly"]]

        query_mat <- x[[cell_type]][["query_mat_subset"]]
        query_anomaly <- x[[cell_type]][["query_anomaly"]]

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
        plot_mat[is.na(plot_mat)] <- 0
        plot_mat[plot_mat > 2] <- 2
        plot_mat[plot_mat < -2] <- -2

        anomaly_labels <- ifelse(anomaly, "Anomalous", "Non-Anomalous")

        annotation_col <- data.frame(
            Dataset = factor(dataset, levels = c("Query", "Reference")),
            Status = factor(anomaly_labels, levels = c("Anomalous", "Non-Anomalous"))
        )
        rownames(annotation_col) <- cell_names

        # Order columns clearly
        order_idx <- order(annotation_col$Dataset, annotation_col$Status)
        plot_mat <- plot_mat[, order_idx, drop = FALSE]
        annotation_col <- annotation_col[order_idx, , drop = FALSE]

        dataset_colors <- c("Query" = "#B565D8", "Reference" = "#5A9BD8")
        anomaly_colors <- c("Non-Anomalous" = "#9E9E9E", "Anomalous" = "#D2314C")

        zscore_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

        top_ha <- ComplexHeatmap::HeatmapAnnotation(
            Status = annotation_col$Status,
            Dataset = annotation_col$Dataset,
            col = list(Status = anomaly_colors, Dataset = dataset_colors),
            show_legend = TRUE,
            show_annotation_name = TRUE,
            annotation_name_side = "left"
        )

        ht <- ComplexHeatmap::Heatmap(
            matrix = plot_mat,
            name = "Z-Score",
            col = zscore_col_fun,
            heatmap_legend_param = list(
                title = "Z-Score",
                at = c(-2, -1, 0, 1, 2),
                labels = c("-2", "-1", "0", "1", "2")
            ),
            show_row_names = nrow(plot_mat) <= 100,
            row_names_gp = if (requireNamespace("grid", quietly = TRUE)) grid::gpar(fontsize = 8) else NULL,
            cluster_rows = TRUE,
            show_column_names = FALSE,
            cluster_columns = FALSE,
            column_order = colnames(plot_mat),
            top_annotation = top_ha,
            border = TRUE,
            column_title = paste0("Local PCA HVG Expression Heatmap: ", cell_type),
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

    # __________________________________________________________
    # PATH B: DISTRIBUTION PLOTS LOGIC (Violin, Boxplot, Ridge)
    # __________________________________________________________

    # Extract threshold and variance info
    threshold <- x[[cell_type]][["applied_threshold"]]
    var_explained <- x[[cell_type]][["var_explained"]]
    total_var <- sum(var_explained)

    # Gather data into a structured data frame
    df_list <- list()

    if (data_type %in% c("reference", "both")) {
        ref_err <- x[[cell_type]][["reference_reconstruction_errors"]]
        ref_anom <- x[[cell_type]][["reference_anomaly"]]

        if (length(ref_err) > 0) {
            df_list[["Reference"]] <- data.frame(
                Error = ref_err,
                Dataset = "Reference",
                Anomaly = ref_anom,
                stringsAsFactors = FALSE
            )
        }
    }

    if (data_type %in% c("query", "both")) {
        if (is.null(x[[cell_type]][["query_reconstruction_errors"]])) {
            if (data_type == "query") {
                stop("There is no query data available in the object to plot.")
            }
        } else {
            query_err <- x[[cell_type]][["query_reconstruction_errors"]]
            query_anom <- x[[cell_type]][["query_anomaly"]]

            if (length(query_err) > 0) {
                df_list[["Query"]] <- data.frame(
                    Error = query_err,
                    Dataset = "Query",
                    Anomaly = query_anom,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(df_list) == 0) {
        stop("No data available to plot for the specified parameters.")
    }

    # Combine and format plotting data
    plot_df <- do.call(rbind, df_list)
    rownames(plot_df) <- NULL

    plot_df[["Dataset"]] <- factor(plot_df[["Dataset"]], levels = c("Reference", "Query"))

    # Updated to Non-Anomalous to match package convention
    plot_df[["Anomaly"]] <- factor(ifelse(plot_df[["Anomaly"]], "Anomalous", "Non-Anomalous"),
                                   levels = c("Non-Anomalous", "Anomalous"))

    # Standardized color palettes for the package
    dataset_colors <- c("Query" = "#B565D8", "Reference" = "#5A9BD8")
    anomaly_colors <- c("Non-Anomalous" = "#9E9E9E", "Anomalous" = "#D2314C")

    # Subtitle with variance info
    sub_title <- sprintf("Threshold: %.2f | PCs utilized explain %.1f%% of Reference variance",
                         threshold, total_var)

    if (plot_type == "violin") {

        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["Dataset"]], y = .data[["Error"]])) +
            ggplot2::geom_violin(ggplot2::aes(fill = .data[["Dataset"]]),
                                 alpha = 0.4, color = "grey50", trim = TRUE, scale = "width") +
            ggplot2::geom_jitter(ggplot2::aes(color = .data[["Anomaly"]]),
                                 width = 0.25, size = 1.5, alpha = 0.7) +
            ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "#D2314C", linewidth = 0.8) +
            ggplot2::scale_fill_manual(values = dataset_colors, guide = "none") +
            ggplot2::scale_color_manual(values = anomaly_colors, name = "Status") +
            ggplot2::labs(
                title = paste0("PCA Reconstruction Error: ", cell_type),
                subtitle = sub_title,
                y = "Reconstruction Error (SSE)",
                x = NULL
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(color = "grey70", linewidth = 0.5),
                axis.text.x = ggplot2::element_text(face = "bold", size = 11),
                axis.title.y = ggplot2::element_text(face = "bold", size = 10),
                legend.position = "right",
                legend.title = ggplot2::element_text(face = "bold", size = 10)
            )

    } else if (plot_type == "boxplot") {

        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["Dataset"]], y = .data[["Error"]])) +
            ggplot2::geom_boxplot(ggplot2::aes(fill = .data[["Dataset"]]),
                                  alpha = 0.5, outlier.shape = NA, width = 0.6) +
            ggplot2::geom_jitter(ggplot2::aes(color = .data[["Anomaly"]]),
                                 width = 0.2, size = 1.5, alpha = 0.6) +
            ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "#D2314C", linewidth = 0.8) +
            ggplot2::scale_fill_manual(values = dataset_colors, guide = "none") +
            ggplot2::scale_color_manual(values = anomaly_colors, name = "Status") +
            ggplot2::labs(
                title = paste0("PCA Reconstruction Error: ", cell_type),
                subtitle = sub_title,
                y = "Reconstruction Error (SSE)",
                x = NULL
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(color = "grey70", linewidth = 0.5),
                axis.text.x = ggplot2::element_text(face = "bold", size = 11),
                axis.title.y = ggplot2::element_text(face = "bold", size = 10),
                legend.position = "right",
                legend.title = ggplot2::element_text(face = "bold", size = 10)
            )

    } else if (plot_type == "ridge") {

        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["Error"]], y = .data[["Dataset"]], fill = .data[["Dataset"]])) +
            ggridges::geom_density_ridges(alpha = 0.6, scale = 1.8, color = "grey30", linewidth = 0.5, rel_min_height = 0.01) +
            ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "#D2314C", linewidth = 0.8) +
            ggplot2::scale_fill_manual(values = dataset_colors, guide = "none") +
            ggplot2::labs(
                title = paste0("PCA Reconstruction Error: ", cell_type),
                subtitle = sub_title,
                x = "Reconstruction Error (SSE)",
                y = NULL
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(color = "grey70", linewidth = 0.5),
                axis.text.y = ggplot2::element_text(face = "bold", size = 11, vjust = 0),
                axis.title.x = ggplot2::element_text(face = "bold", size = 10)
            )
    }

    if (draw_plot && plot_type != "heatmap") {
        print(p)
        return(invisible(p))
    } else {
        return(p)
    }
}
