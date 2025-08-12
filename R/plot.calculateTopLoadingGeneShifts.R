#' @title Plot Top Loading Gene Expression Shifts
#'
#' @description
#' This function creates visualizations showing expression distributions for top loading genes
#' that exhibit distributional differences between query and reference datasets. Can display
#' results as elegant complex heatmaps or as information-rich summary boxplots.
#'
#' @details
#' This function visualizes the results from \code{calculateTopLoadingGeneShifts}.
#' The "heatmap" option displays a hierarchically clustered set of genes.
#' The "boxplot" option creates a two-panel plot using `ggplot2`: the left panel shows
#' horizontal expression boxplots for up to 5 PCs, while the right panel displays their
#' corresponding PC loadings and adjusted p-values.
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject}.
#' @param cell_type A character string specifying the cell type to plot (must be exactly one).
#' @param pc_subset A numeric vector specifying which principal components to plot. Default is 1:3.
#' @param plot_type A character string specifying visualization type. Either "heatmap" or "boxplot".
#'                  Default is "heatmap".
#' @param plot_by A character string specifying gene selection method when `n_genes` is not NULL.
#'                Either "top_loading" or "p_adjusted". Default is "p_adjusted".
#' @param n_genes Number of top genes to show per PC. Can be NULL if `significance_threshold` is set.
#'                Default is 10.
#' @param significance_threshold If not NULL, a numeric value between 0 and 1. Used for gene
#'   selection or annotation. Default is 0.05.
#' @param ... Additional arguments passed to \code{\link[ComplexHeatmap]{draw}} or not used for boxplot.
#'
#' @return A plot object.
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateTopLoadingGeneShifts}}
#'
#' @rdname calculateTopLoadingGeneShifts
#'
#' @importFrom stats wilcox.test var p.adjust na.omit setNames aggregate reshape
#'
# Function to visualize genes with top loadings
plot.calculateTopLoadingGeneShiftsObject <- function(x,
                                                     cell_type,
                                                     pc_subset = 1:3,
                                                     plot_type = c("heatmap", "boxplot"),
                                                     plot_by = c("p_adjusted", "top_loading"),
                                                     n_genes = 10,
                                                     significance_threshold = 0.05,
                                                     ...) {

    # --- Input Validation ---
    if (missing(cell_type) || length(cell_type) != 1) {
        stop("cell_type must be specified and be a single character string.")
    }

    # Match arguments first
    plot_type <- match.arg(plot_type)
    plot_by <- match.arg(plot_by)

    # Validate plot-type specific requirements
    if (plot_type == "boxplot" && is.null(n_genes)) {
        stop("For 'boxplot' plot type, 'n_genes' must be specified.")
    }
    if (plot_type == "heatmap" &&
        is.null(n_genes) &&
        is.null(significance_threshold)) {
        stop("For 'heatmap' plot type, at least one of 'n_genes' or 'significance_threshold' must be specified.")
    }
    if (!is.null(n_genes) &&
        (!is.numeric(n_genes) ||
         length(n_genes) != 1 ||
         n_genes <= 0)) {
        stop("If not NULL, 'n_genes' must be a positive integer.")
    }
    if (!is.null(significance_threshold) &&
        (!is.numeric(significance_threshold) ||
         length(significance_threshold) != 1 ||
         significance_threshold < 0 ||
         significance_threshold > 1)) {
        stop("If not NULL, 'significance_threshold' must be a numeric value between 0 and 1.")
    }

    # Check for package dependencies
    if (plot_type == "heatmap") {
        if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
            stop("Packages 'ComplexHeatmap' and 'circlize' are required. Please install them.",
                 call. = FALSE)
        }
    }

    # Check object structure
    required_elements <- c("expression_data", "cell_metadata", "percent_var")
    if (!all(required_elements %in% names(x))) {
        stop("Input object 'x' is missing required elements.", call. = FALSE)
    }

    # Get available PCs from the object
    pc_names <- paste0("PC", pc_subset)
    available_pcs <- intersect(pc_names, names(x))
    if (length(available_pcs) == 0) {
        stop("None of the requested PCs were found in the results object.")
    }
    available_pcs <- pc_names[pc_names %in% available_pcs]

    # Limit PCs for boxplot
    if (plot_type == "boxplot" && length(available_pcs) > 5) {
        warning("Boxplot type can display at most 5 PCs. Using the first 5.")
        available_pcs <- available_pcs[1:5]
    }

    # Check cell type existence
    if (!cell_type %in% x[["cell_metadata"]][["cell_type"]]) {
        stop(paste("Cell type '", cell_type, "' not found in results.", sep = ""))
    }

    # --- Route to Plotting Functions ---
    if (plot_type == "boxplot") {
        return(plotBoxplot(x, cell_type, available_pcs,
                           plot_by, n_genes, significance_threshold))
    } else {
        ht <- plotHeatmap(x, cell_type, available_pcs,
                          plot_by, n_genes, significance_threshold)
        if (is.null(ht)) {
            message("No data available to generate a heatmap for the specified parameters.")
            return(invisible(NULL))
        }
        # Draw the heatmap
        return(
            ComplexHeatmap::draw(ht,
                                 heatmap_legend_side = "right",
                                 annotation_legend_side = "right",
                                 merge_legends = TRUE,
                                 ...)
        )
    }
}

#' @title Plot Heatmaps for Top Loading Gene Shifts (Simplified Single Heatmap)
#'
#' @description
#' This internal helper function creates a single, hierarchically clustered heatmap
#' displaying expression data for top loading genes from principal component analysis.
#' The function handles gene selection, data preprocessing, and visualization formatting.
#'
#' @details
#' This function generates a ComplexHeatmap visualization showing scaled gene expression
#' data across cells. Genes are selected based on either their loading values or
#' statistical significance. The function automatically handles data filtering,
#' scaling, and visual formatting including color schemes and annotations.
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject} containing
#'   expression data and analysis results.
#' @param cell_type A character string specifying the cell type to visualize.
#' @param available_pcs A character vector of principal components to include in analysis.
#' @param plot_by A character string indicating gene selection criterion ("top_loading" or "p_adjusted").
#' @param n_genes An integer specifying the number of top genes to display per PC. Can be NULL.
#' @param significance_threshold A numeric value between 0 and 1 for significance filtering. Can be NULL.
#'
#' @keywords internal
#'
#' @return A ComplexHeatmap object ready for plotting, or NULL if no genes meet selection criteria.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Helper heatmap function
plotHeatmap <- function(x, cell_type, available_pcs, plot_by,
                        n_genes, significance_threshold) {

    # Initialize gene collection
    all_genes_to_plot <- c()

    # Gene selection based on n_genes criterion
    if (!is.null(n_genes)) {
        gene_list <- lapply(available_pcs, function(pc_name) {
            pc_results <- x[[pc_name]]
            if (is.null(pc_results) || nrow(pc_results) == 0) {
                return(NULL)
            }

            pc_cell_data <- pc_results[pc_results[["cell_type"]] == cell_type, ]
            if (nrow(pc_cell_data) == 0) {
                return(NULL)
            }

            # Order genes by selection criterion
            if (plot_by == "top_loading") {
                gene_order <- order(abs(pc_cell_data[["loading"]]),
                                    decreasing = TRUE)
            } else {
                gene_order <- order(pc_cell_data[["p_adjusted"]])
            }

            selected_indices <- gene_order[1:min(n_genes, nrow(pc_cell_data))]
            return(pc_cell_data[selected_indices, "gene"])
        })
        all_genes_to_plot <- unique(unlist(gene_list))

    } else if (!is.null(significance_threshold)) {
        # Gene selection based on significance threshold
        sig_gene_list <- lapply(available_pcs, function(pc_name) {
            pc_results <- x[[pc_name]]
            if (is.null(pc_results) || nrow(pc_results) == 0) {
                return(NULL)
            }

            sig_df <- pc_results[pc_results[["cell_type"]] == cell_type &
                                     pc_results[["p_adjusted"]] < significance_threshold, ]
            return(sig_df[["gene"]])
        })
        all_genes_to_plot <- unique(unlist(sig_gene_list))
    }

    # Check if any genes were selected
    if (length(all_genes_to_plot) == 0) {
        return(NULL)  # Return NULL instead of warning
    }

    # Identify significant genes for annotation
    significant_genes <- c()
    if (!is.null(significance_threshold)) {
        sig_gene_list <- lapply(available_pcs, function(pc_name) {
            pc_results <- x[[pc_name]]
            if (is.null(pc_results) || nrow(pc_results) == 0) {
                return(NULL)
            }

            sig_df <- pc_results[pc_results[["cell_type"]] == cell_type &
                                     pc_results[["p_adjusted"]] < significance_threshold, ]
            return(sig_df[["gene"]])
        })
        significant_genes <- unique(unlist(sig_gene_list))
    }

    # Prepare cell subset and ordering
    cell_subset <- x[["cell_metadata"]][x[["cell_metadata"]][["cell_type"]] == cell_type, ]
    cell_subset <- cell_subset[order(cell_subset[["dataset"]]), ]

    # Extract and filter expression matrix
    expr_matrix_unfiltered <- as.matrix(
        x[["expression_data"]][all_genes_to_plot,
                               cell_subset[["cell_id"]],
                               drop = FALSE]
    )

    # Remove zero-variance genes to prevent scaling issues
    row_variances <- apply(expr_matrix_unfiltered, 1, stats::var)
    genes_to_keep <- row_variances > .Machine$double.eps

    if (sum(genes_to_keep) == 0) {
        return(NULL)  # Return NULL instead of warning
    }

    # Create final expression matrix and scale
    expr_matrix <- expr_matrix_unfiltered[genes_to_keep, , drop = FALSE]
    scaled_matrix <- t(scale(t(expr_matrix)))

    # Clamp z-scores to the range [-2, 2] for consistent visualization
    scaled_matrix[scaled_matrix < -2] <- -2
    scaled_matrix[scaled_matrix > 2] <- 2

    # Define color schemes
    dataset_colors <- c("Query" = "#E41A1C", "Reference" = "#377EB8")

    # Create custom color mapping function
    .createColorMapping <- function(breaks, colors) {
        function(x) {
            # Handle NA values
            if (any(is.na(x))) {
                return(rep("grey", length(x)))
            }

            # Clamp values to break range
            x[x < min(breaks)] <- min(breaks)
            x[x > max(breaks)] <- max(breaks)

            # Create color ramp and normalize values
            col_ramp <- grDevices::colorRamp(colors)
            x_norm <- (x - min(breaks)) / (max(breaks) - min(breaks))

            # Convert to hex colors
            rgb_vals <- col_ramp(x_norm)
            grDevices::rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue = 255)
        }
    }
    zscore_col_fun <- .createColorMapping(c(-2, 0, 2),
                                          c("#313695", "white", "#A50026"))

    # Create top annotation
    top_ha <- ComplexHeatmap::HeatmapAnnotation(
        Dataset = cell_subset[["dataset"]],
        col = list(Dataset = dataset_colors),
        show_legend = TRUE,
        show_annotation_name = FALSE,
        annotation_legend_param = list(title = "Dataset")
    )

    # Prepare gene labels with significance indicators
    final_genes <- rownames(expr_matrix)
    row_labels <- paste0(final_genes,
                         ifelse(final_genes %in% significant_genes, "*", ""))

    # Smart dependency checking for grid package
    if (requireNamespace("grid", quietly = TRUE)) {
        # Use smaller font size for better readability when grid is available
        row_gp <- grid::gpar(fontsize = 8)
        ht <- ComplexHeatmap::Heatmap(
            matrix = scaled_matrix,
            name = "Z-Score",
            col = zscore_col_fun,
            heatmap_legend_param = list(
                title = "Z-Score",
                at = c(-2, -1, 0, 1, 2),
                labels = c("-2", "-1", "0", "1", "2")
            ),
            show_row_names = TRUE,
            row_labels = row_labels,
            row_names_gp = row_gp,
            cluster_rows = TRUE,
            show_column_names = FALSE,
            cluster_columns = FALSE,
            column_order = cell_subset[["cell_id"]],
            top_annotation = top_ha,
            border = TRUE,
            use_raster = TRUE,
            raster_quality = 2
        )
    } else {
        # Use default font size and inform user about grid package
        message("Note: For optimal gene name readability, consider installing the 'grid' package: install.packages('grid')")
        ht <- ComplexHeatmap::Heatmap(
            matrix = scaled_matrix,
            name = "Z-Score",
            col = zscore_col_fun,
            heatmap_legend_param = list(
                title = "Z-Score",
                at = c(-2, -1, 0, 1, 2),
                labels = c("-2", "-1", "0", "1", "2")
            ),
            show_row_names = TRUE,
            row_labels = row_labels,
            cluster_rows = TRUE,
            show_column_names = FALSE,
            cluster_columns = FALSE,
            column_order = cell_subset[["cell_id"]],
            top_annotation = top_ha,
            border = TRUE,
            use_raster = TRUE,
            raster_quality = 2
        )
    }

    return(ht)
}

#' @title Plot Boxplots for Top Loading Gene Shifts (Two-Panel Summary using ggplot2)
#'
#' @description
#' This internal helper function creates a comprehensive two-panel summary plot
#' displaying gene expression distributions and principal component loadings.
#' The visualization uses ggplot2 faceting to create side-by-side panels.
#'
#' @details
#' This function generates a dual-panel ggplot2 visualization where the left panel
#' shows horizontal boxplots of gene expression distributions comparing Reference
#' and Query datasets, while the right panel displays PC loading values as points
#' with adjusted p-values. Gene selection is based on the union of top genes
#' across specified principal components.
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject} containing
#'   expression data and analysis results.
#' @param cell_type A character string specifying the cell type to visualize.
#' @param available_pcs A character vector of principal components to include in analysis.
#' @param plot_by A character string indicating gene selection criterion ("top_loading" or "p_adjusted").
#' @param n_genes An integer specifying the number of top genes to display per PC.
#' @param significance_threshold A numeric value between 0 and 1 for significance annotation.
#'
#' @keywords internal
#'
#' @return A ggplot2 object ready for display, or NULL if no genes meet selection criteria.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Helper boxplot function - unchanged but with line breaks
plotBoxplot <- function(x, cell_type, available_pcs, plot_by,
                        n_genes, significance_threshold) {

    # Validate required parameters
    if (is.null(n_genes)) {
        stop("'plot_type = \"boxplot\"' requires 'n_genes' to be set.")
    }

    # Gene selection across principal components
    gene_df_list <- lapply(available_pcs, function(pc_name) {
        pc_results <- x[[pc_name]]
        if (is.null(pc_results) || nrow(pc_results) == 0) {
            return(NULL)
        }

        pc_cell_data <- pc_results[pc_results[["cell_type"]] == cell_type, ]
        if (nrow(pc_cell_data) == 0) {
            return(NULL)
        }

        # Order genes by selection criterion
        if (plot_by == "top_loading") {
            gene_order <- order(abs(pc_cell_data[["loading"]]),
                                decreasing = TRUE)
        } else {
            gene_order <- order(pc_cell_data[["p_adjusted"]])
        }

        selected_indices <- gene_order[1:min(n_genes, nrow(pc_cell_data))]
        return(pc_cell_data[selected_indices, ])
    })

    # Combine results and remove duplicates
    final_gene_data <- do.call(rbind, gene_df_list)
    final_gene_data <- final_gene_data[!duplicated(final_gene_data[["gene"]]), ]

    if (nrow(final_gene_data) == 0) {
        warning("No genes met the selection criteria for the boxplot.")
        return(NULL)
    }

    # Order final gene set by selection criterion
    if (plot_by == "top_loading") {
        final_gene_data <- final_gene_data[
            order(abs(final_gene_data[["loading"]]), decreasing = TRUE),
        ]
    } else {
        final_gene_data <- final_gene_data[
            order(final_gene_data[["p_adjusted"]]),
        ]
    }

    # Prepare expression data for plotting
    cell_subset_meta <- x[["cell_metadata"]][
        x[["cell_metadata"]][["cell_type"]] == cell_type,
    ]
    expr_subset <- x[["expression_data"]][
        final_gene_data[["gene"]],
        cell_subset_meta[["cell_id"]],
        drop = FALSE
    ]

    # Create expression data frame
    expr_plot_list <- lapply(final_gene_data[["gene"]], function(g) {
        data.frame(
            gene = g,
            value = as.numeric(expr_subset[g, ]),
            dataset = cell_subset_meta[["dataset"]],
            metric = "Expression",
            pc_group = NA_character_,
            stringsAsFactors = FALSE
        )
    })
    expr_plot_df <- do.call(rbind, expr_plot_list)

    # Prepare loading data for plotting
    loading_data <- x[["gene_metadata"]][
        x[["gene_metadata"]][["gene"]] %in% final_gene_data[["gene"]] &
            paste0("PC", x[["gene_metadata"]][["pc"]]) %in% available_pcs,
    ]

    loading_plot_df <- data.frame(
        gene = loading_data[["gene"]],
        value = loading_data[["loading"]],
        dataset = factor("Reference", levels = c("Reference", "Query")),
        metric = "PC Loading",
        pc_group = paste0("PC", loading_data[["pc"]]),
        stringsAsFactors = FALSE
    )

    # Combine plotting data
    full_plot_df <- rbind(expr_plot_df, loading_plot_df)

    # Prepare gene labels with significance indicators
    is_significant <- if (!is.null(significance_threshold)) {
        final_gene_data[["p_adjusted"]] < significance_threshold
    } else {
        FALSE
    }
    gene_labels <- paste0(final_gene_data[["gene"]],
                          ifelse(is_significant, "*", ""))
    full_plot_df[["gene"]] <- factor(full_plot_df[["gene"]],
                                     levels = rev(final_gene_data[["gene"]]))

    # Format p-value labels
    p_value_data <- final_gene_data[, c("gene", "p_adjusted")]
    p_value_labels <- paste0("p = ",
                             sprintf("%.2g", p_value_data[["p_adjusted"]]))
    p_value_labels[p_value_data[["p_adjusted"]] < 0.001] <- "p < .001"
    names(p_value_labels) <- p_value_data[["gene"]]

    # Create PC legend labels with variance explained
    pc_legend_labels <- paste0(available_pcs, " (",
                               sprintf("%.1f%%", x[["percent_var"]][available_pcs]),
                               ")")
    names(pc_legend_labels) <- available_pcs

    # Define color and shape palettes
    dataset_colors <- c("Reference" = "#377EB8", "Query" = "#E41A1C")
    num_pcs <- length(available_pcs)
    pc_color_palette <- grDevices::hcl.colors(num_pcs, palette = "Dark 3")
    pc_shape_palette <- c(16, 17, 15, 18, 8)

    pc_colors <- pc_color_palette[1:num_pcs]
    pc_shapes <- pc_shape_palette[1:num_pcs]
    names(pc_colors) <- available_pcs
    names(pc_shapes) <- available_pcs

    # Create facet labeller
    facet_labeller <- ggplot2::labeller(
        metric = c("Expression" = "Log-Normalized Expression",
                   "PC Loading" = "PC Loading")
    )

    # Build the ggplot object
    p <- ggplot2::ggplot(
        mapping = ggplot2::aes(y = .data[["gene"]], x = .data[["value"]])
    ) +

        # Expression panel: horizontal boxplots
        ggplot2::geom_boxplot(
            data = ~subset(full_plot_df, metric == "Expression"),
            ggplot2::aes(fill = .data[["dataset"]]),
            alpha = 0.8,
            outlier.size = 0.5
        ) +

        # Loading panel: reference line and points
        ggplot2::geom_vline(xintercept = 0,
                            linetype = "dashed",
                            color = "grey50") +
        ggplot2::geom_point(
            data = ~subset(full_plot_df, metric == "PC Loading"),
            ggplot2::aes(color = .data[["pc_group"]],
                         shape = .data[["pc_group"]]),
            size = 3,
            alpha = 0.8
        ) +

        # Panel configuration
        ggplot2::facet_grid(~ metric,
                            scales = "free_x",
                            switch = "x",
                            labeller = facet_labeller) +

        # Scale definitions
        ggplot2::scale_fill_manual(name = "Dataset", values = dataset_colors) +
        ggplot2::scale_color_manual(name = "PC",
                                    values = pc_colors,
                                    labels = pc_legend_labels) +
        ggplot2::scale_shape_manual(name = "PC",
                                    values = pc_shapes,
                                    labels = pc_legend_labels) +
        ggplot2::scale_y_discrete(labels = rev(gene_labels), name = NULL) +
        ggplot2::scale_x_continuous(
            sec.axis = ggplot2::dup_axis(
                name = "Adj. P-value",
                breaks = NULL,
                labels = rev(p_value_labels)
            ),
            name = NULL,
            expand = ggplot2::expansion(mult = c(0.05, 0.15))
        ) +

        # Labels and guides
        ggplot2::labs(title = NULL, subtitle = NULL) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(order = 1),
            color = ggplot2::guide_legend(order = 2),
            shape = ggplot2::guide_legend(order = 2)
        ) +

        # Theme customization
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.y.left = ggplot2::element_text(face = "bold"),
            legend.position = "right",
            legend.box = "vertical",
            strip.background = ggplot2::element_blank(),
            strip.placement = "outside",
            strip.text.x = ggplot2::element_text(face = "bold", size = 10)
        )

    return(p)
}
