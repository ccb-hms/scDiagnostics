#' @title Plot Top Loading Gene Expression Shifts
#'
#' @description
#' This function creates visualizations showing expression distributions for top loading genes
#' that exhibit distributional differences between query and reference datasets. Can display
#' results as elegant complex heatmaps, information-rich summary boxplots, or pseudo-bulk fold
#' change barplots. Optionally displays anomaly status when available.
#'
#' @details
#' This function visualizes the results from \code{calculateTopLoadingGeneShifts}.
#' The "heatmap" option displays a hierarchically clustered set of genes.
#' The "boxplot" option creates a two-panel plot using `ggplot2`: the left panel shows
#' horizontal expression boxplots for up to 5 PCs, while the right panel displays their
#' corresponding PC loadings and adjusted p-values.
#' The "barplot" option creates horizontal barplots showing log2 fold changes between
#' pseudo-bulk expression profiles (query vs reference), with genes ordered identically
#' to the heatmap clustering. Bars show comparisons for query non-anomaly (green),
#' optionally all query cells (yellow), and query anomaly cells (red) versus reference.
#' When anomaly detection results are available and \code{show_anomalies} is TRUE,
#' additional annotation bars or visual cues highlight anomalous cells.
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject}.
#' @param cell_type A character string specifying the cell type to plot (must be exactly one).
#' @param pc_subset A numeric vector specifying which principal components to plot. Default is 1:3.
#' @param plot_type A character string specifying visualization type. Either "heatmap", "barplot", or "boxplot".
#'                  Default is "heatmap".
#' @param plot_by A character string specifying gene selection method when `n_genes` is not NULL.
#'                Either "top_loading" or "p_adjusted". Default is "p_adjusted".
#' @param n_genes Number of top genes to show per PC. Can be NULL if `significance_threshold` is set.
#'                Default is 10.
#' @param significance_threshold If not NULL, a numeric value between 0 and 1. Used for gene
#'   selection or annotation. Default is 0.05.
#' @param show_anomalies Logical indicating whether to display anomaly status annotations.
#'                      Default is FALSE. Requires anomaly results to be present in the object.
#' @param show_all_query Logical indicating whether to show the yellow bar representing all
#'                       query cells vs reference in barplot visualization. Only applies when
#'                       \code{plot_type = "barplot"} and anomaly data is available. Default is TRUE.
#' @param pseudo_bulk Logical indicating whether to create pseudo-bulk profiles instead of
#'                    showing individual cells. When TRUE, expression values are averaged within groups
#'                    (dataset and optionally anomaly status). Not compatible with boxplot visualization.
#'                    Required for barplot visualization. Default is FALSE.
#' @param cluster_cols Logical indicating whether to cluster columns in the heatmap when
#'                    `pseudo_bulk = TRUE`. When TRUE, columns (pseudo-bulk profiles) will be
#'                    hierarchically clustered. When FALSE, columns maintain their original ordering
#'                    (Query groups followed by Reference groups). Only applicable when
#'                    `pseudo_bulk = TRUE` and `plot_type = "heatmap"`. Default is FALSE.
#' @param draw_plot Logical indicating whether to draw the plot immediately (TRUE) or return
#'                  the undrawn plot object (FALSE). For heatmaps, FALSE returns a ComplexHeatmap
#'                  object that can be further customized before drawing. Default is TRUE.
#' @param show_all_query Logical indicating whether to show the yellow bar for all query vs reference
#'                       comparison. Default is TRUE. When FALSE, only green and red bars are shown.
#' @param max_cells_ref Maximum number of reference cells to include in the plot. If NULL,
#' all available reference cells are plotted. Default is NULL.
#' @param max_cells_query Maximum number of query cells to include in the plot. If NULL,
#' all available query cells are plotted. Default is NULL.
#' @param ... Additional arguments passed to \code{\link[ComplexHeatmap]{draw}} or not used for other plot types.
#'
#' @return A plot object. For heatmaps when \code{draw_plot = FALSE}, returns a ComplexHeatmap object.
#' For boxplots and barplots, returns a ggplot2 object.
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
                                                     plot_type = c("heatmap", "barplot", "boxplot"),
                                                     plot_by = c("p_adjusted", "top_loading"),
                                                     n_genes = 10,
                                                     significance_threshold = 0.05,
                                                     show_anomalies = FALSE,
                                                     pseudo_bulk = FALSE,
                                                     cluster_cols = FALSE,
                                                     draw_plot = TRUE,
                                                     show_all_query = TRUE,
                                                     max_cells_ref = NULL,
                                                     max_cells_query = NULL,
                                                     ...) {

    #Input Validation
    if (missing(cell_type) || length(cell_type) != 1) {
        stop("cell_type must be specified and be a single character string.")
    }

    # Match arguments first
    plot_type <- match.arg(plot_type)
    plot_by <- match.arg(plot_by)

    # Validate show_anomalies parameter
    if (!is.logical(show_anomalies) || length(show_anomalies) != 1) {
        stop("show_anomalies must be a logical value.")
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

    # Add after existing validation
    if (!is.logical(pseudo_bulk) || length(pseudo_bulk) != 1) {
        stop("pseudo_bulk must be a logical value.")
    }

    if (!is.logical(cluster_cols) || length(cluster_cols) != 1) {
        stop("cluster_cols must be a logical value.")
    }

    # Restrict boxplot with pseudo_bulk
    if (pseudo_bulk && plot_type == "boxplot") {
        stop("Boxplot visualization is not compatible with pseudo_bulk = TRUE. ",
             "Please use plot_type = 'heatmap' or set pseudo_bulk = FALSE.")
    }

    # Require pseudo_bulk for barplot
    if (plot_type == "barplot" && !pseudo_bulk) {
        stop("Barplot visualization requires pseudo_bulk = TRUE. ",
             "Please set pseudo_bulk = TRUE for barplot.")
    }

    # Restrict cluster_cols without pseudo_bulk
    if (cluster_cols && !pseudo_bulk) {
        stop("cluster_cols = TRUE requires pseudo_bulk = TRUE.")
    }

    # Check if anomaly data is available when requested
    if (show_anomalies) {
        if (!"anomaly_status" %in% names(x[["cell_metadata"]])) {
            stop("Anomaly visualization requested but anomaly data not found in object. ",
                 "Please re-run calculateTopLoadingGeneShifts() with detect_anomalies = TRUE.")
        }
    }

    # Validate plot-type specific requirements
    if (plot_type %in% c("boxplot", "barplot") && is.null(n_genes)) {
        stop("For '", plot_type, "' plot type, 'n_genes' must be specified.")
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

    # Downsample cells if max_cells parameters are specified
    cell_metadata_filtered <- x[["cell_metadata"]]
    cell_subset <- cell_metadata_filtered[cell_metadata_filtered[["cell_type"]] == cell_type, ]

    # Separate reference and query data for potential downsampling
    ref_cells <- cell_subset[cell_subset[["dataset"]] == "Reference", ]
    query_cells <- cell_subset[cell_subset[["dataset"]] == "Query", ]

    # Downsample reference cells if specified
    if (!is.null(max_cells_ref) && nrow(ref_cells) > max_cells_ref) {
        sampled_indices <- sample(nrow(ref_cells), max_cells_ref)
        ref_cells <- ref_cells[sampled_indices, ]
    }

    # Downsample query cells if specified
    if (!is.null(max_cells_query) && nrow(query_cells) > max_cells_query) {
        sampled_indices <- sample(nrow(query_cells), max_cells_query)
        query_cells <- query_cells[sampled_indices, ]
    }

    # Combine the potentially downsampled cells back together
    cell_subset_final <- rbind(ref_cells, query_cells)

    # Update the object with downsampled data for plotting
    x_plotting <- x
    x_plotting[["cell_metadata"]] <- cell_subset_final
    x_plotting[["expression_data"]] <- x[["expression_data"]][, cell_subset_final[["cell_id"]],
                                                              drop = FALSE]

    # Route to Plotting Functions
    if (plot_type == "heatmap") {
        ht <- plotHeatmap(x_plotting, cell_type,
                          available_pcs, plot_by, n_genes,
                          significance_threshold, show_anomalies,
                          pseudo_bulk, cluster_cols)

        if (is.null(ht)) {
            message("No data available to generate a heatmap for the specified parameters.")
            return(invisible(NULL))
        }

        # Return drawn or undrawn heatmap based on draw_plot parameter
        if (draw_plot) {
            # Draw the heatmap and return the drawn object
            return(
                ComplexHeatmap::draw(ht,
                                     heatmap_legend_side = "right",
                                     annotation_legend_side = "right",
                                     merge_legends = TRUE,
                                     ...)
            )
        } else {
            # Return the undrawn Heatmap object for further customization
            return(ht)
        }
    } else if (plot_type == "barplot") {
        return(plotBarplot(x_plotting, cell_type,
                           available_pcs, plot_by,
                           n_genes, significance_threshold,
                           show_anomalies, show_all_query))
    } else {
        return(plotBoxplot(x_plotting, cell_type,
                           available_pcs, plot_by,
                           n_genes, significance_threshold,
                           show_anomalies))
    }
}

#' @title Plot Heatmaps for Top Loading Gene Shifts (Simplified Single Heatmap)
#'
#' @description
#' This internal helper function creates a single, hierarchically clustered heatmap
#' displaying expression data for top loading genes from principal component analysis.
#' The function handles gene selection, data preprocessing, and visualization formatting.
#' Optionally includes anomaly status annotations with proper cell ordering.
#'
#' @details
#' This function generates a ComplexHeatmap visualization showing scaled gene expression
#' data across cells. Genes are selected based on either their loading values or
#' statistical significance. The function automatically handles data filtering,
#' scaling, and visual formatting including color schemes and annotations.
#'
#' Cell ordering (left to right):
#' 1. Query Anomalous cells (leftmost)
#' 2. Query Normal cells
#' 3. Reference Anomalous cells
#' 4. Reference Normal cells (rightmost)
#'
#' This ensures clear dataset separation while grouping anomalous cells together
#' within each dataset for easy visual identification.
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject} containing
#'          expression data and analysis results.
#' @param cell_type A character string specifying the cell type to visualize.
#' @param available_pcs A character vector of principal components to include in analysis.
#' @param plot_by A character string indicating gene selection criterion ("top_loading" or "p_adjusted").
#' @param n_genes An integer specifying the number of top genes to display per PC. Can be NULL.
#' @param significance_threshold A numeric value between 0 and 1 for significance filtering. Can be NULL.
#' @param show_anomalies Logical indicating whether to show anomaly annotations.
#' @param pseudo_bulk Logical indicating whether to create pseudo-bulk profiles instead of
#'                    showing individual cells. When TRUE, expression values are averaged within groups
#'                    (dataset and optionally anomaly status). Not compatible with boxplot visualization.
#'                    Default is FALSE.
#' @param cluster_cols Logical indicating whether to cluster columns in the heatmap when
#'                    `pseudo_bulk = TRUE`. When TRUE, columns (pseudo-bulk profiles) will be
#'                    hierarchically clustered. When FALSE, columns maintain their original ordering
#'                    (Query groups followed by Reference groups). Only applicable when
#'                    `pseudo_bulk = TRUE` and `plot_type = "heatmap"`. Default is FALSE.
#'
#' @keywords internal
#'
#' @return A ComplexHeatmap object ready for plotting, or NULL if no genes meet selection criteria.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Helper heatmap function
plotHeatmap <- function(x, cell_type,
                        available_pcs, plot_by,
                        n_genes, significance_threshold,
                        show_anomalies, pseudo_bulk = FALSE,
                        cluster_cols = FALSE) {

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

    # Prepare cell subset
    cell_subset <- x[["cell_metadata"]][x[["cell_metadata"]][["cell_type"]] == cell_type, ]

    # Handle pseudo-bulk aggregation
    if (pseudo_bulk) {
        # Create grouping categories
        if (show_anomalies && "anomaly_status" %in% names(cell_subset)) {
            # Only show anomalies for query cells, treat all reference as normal
            cell_subset_modified <- cell_subset
            cell_subset_modified[["anomaly_status"]][cell_subset_modified[["dataset"]] == "Reference"] <- "Normal"

            # 3 categories: Reference_Normal, Query_Normal, Query_Anomaly (no Reference_Anomaly)
            cell_subset_modified[["group"]] <- paste(cell_subset_modified[["dataset"]],
                                                     cell_subset_modified[["anomaly_status"]],
                                                     sep = "_")
            group_order <- c("Reference_Normal", "Query_Normal", "Query_Anomaly")

            # Update cell_subset to use modified version
            cell_subset <- cell_subset_modified
            cell_subset[["group"]] <- cell_subset_modified[["group"]]
        } else {
            # 2 categories: Reference, Query
            cell_subset[["group"]] <- cell_subset[["dataset"]]
            group_order <- c("Reference", "Query")
        }

        # Aggregate expression data by groups
        expr_matrix_unfiltered <- as.matrix(
            x[["expression_data"]][all_genes_to_plot, cell_subset[["cell_id"]], drop = FALSE]
        )

        # Create pseudo-bulk profiles
        pseudo_bulk_matrix <- matrix(0, nrow = nrow(expr_matrix_unfiltered), ncol = length(group_order))
        rownames(pseudo_bulk_matrix) <- rownames(expr_matrix_unfiltered)
        colnames(pseudo_bulk_matrix) <- group_order

        for (group in group_order) {
            group_cells <- cell_subset[cell_subset[["group"]] == group, "cell_id"]
            if (length(group_cells) > 0) {
                if (length(group_cells) == 1) {
                    pseudo_bulk_matrix[, group] <- expr_matrix_unfiltered[, group_cells]
                } else {
                    pseudo_bulk_matrix[, group] <- rowMeans(expr_matrix_unfiltered[, group_cells, drop = FALSE])
                }
            }
        }

        # Update expression matrix for plotting
        expr_matrix_unfiltered <- pseudo_bulk_matrix

        # Create new cell_subset for annotations
        cell_subset_pseudo <- data.frame(
            cell_id = group_order,
            group = group_order,
            stringsAsFactors = FALSE
        )

        # Parse group names back to dataset and anomaly status
        if (show_anomalies && "anomaly_status" %in% names(cell_subset)) {
            cell_subset_pseudo[["dataset"]] <- gsub("_.*", "", cell_subset_pseudo[["group"]])
            cell_subset_pseudo[["anomaly_status"]] <- gsub(".*_", "", cell_subset_pseudo[["group"]])

            # Ensure reference cells are marked as Normal (redundant safety check)
            cell_subset_pseudo[["anomaly_status"]][cell_subset_pseudo[["dataset"]] == "Reference"] <- "Normal"
        } else {
            cell_subset_pseudo[["dataset"]] <- cell_subset_pseudo[["group"]]
        }

        cell_subset <- cell_subset_pseudo
    } else {
        # Original ordering logic for individual cells
        if (show_anomalies && "anomaly_status" %in% names(cell_subset)) {
            # Only show anomalies for query cells, treat all reference as normal
            cell_subset[["anomaly_status"]][cell_subset[["dataset"]] == "Reference"] <- "Normal"

            cell_subset <- cell_subset[order(
                factor(cell_subset[["dataset"]], levels = c("Reference", "Query")),
                factor(cell_subset[["anomaly_status"]], levels = c("Normal", "Anomaly"))
            ), ]
        } else {
            cell_subset <- cell_subset[order(
                factor(cell_subset[["dataset"]], levels = c("Reference", "Query"))
            ), ]
        }

        # Extract expression matrix for individual cells
        expr_matrix_unfiltered <- as.matrix(
            x[["expression_data"]][all_genes_to_plot, cell_subset[["cell_id"]], drop = FALSE]
        )
    }

    # Remove zero-variance genes
    row_variances <- apply(expr_matrix_unfiltered, 1, stats::var)
    genes_to_keep <- row_variances > .Machine$double.eps

    if (sum(genes_to_keep) == 0) {
        return(NULL)
    }

    # Create final expression matrix and scale
    expr_matrix <- expr_matrix_unfiltered[genes_to_keep, , drop = FALSE]
    scaled_matrix <- t(scale(t(expr_matrix)))

    # Clamp z-scores to the range [-2, 2] for consistent visualization
    scaled_matrix[scaled_matrix < -2] <- -2
    scaled_matrix[scaled_matrix > 2] <- 2

    # Define color schemes
    dataset_colors <- c("Query" = "#B565D8", "Reference" = "#5A9BD8")
    anomaly_colors <- c("Non-Anomalous" = "#9E9E9E", "Anomalous" = "#D2314C")

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

    # Create top annotation with query-only anomaly status
    annotation_list <- list()
    color_list <- list()
    legend_params <- list()

    if (show_anomalies && "anomaly_status" %in% names(cell_subset)) {
        # Map anomaly status to new labels (reference cells are already set to "Normal")
        anomaly_labels <- ifelse(cell_subset[["anomaly_status"]] == "Anomaly",
                                 "Anomalous", "Non-Anomalous")
        annotation_list[["Status"]] <- anomaly_labels
        color_list[["Status"]] <- anomaly_colors
        legend_params[["Status"]] <- list(title = "Query Anomaly Status")
    }

    # Add dataset annotation below anomaly status
    annotation_list[["Dataset"]] <- cell_subset[["dataset"]]
    color_list[["Dataset"]] <- dataset_colors
    legend_params[["Dataset"]] <- list(title = "Dataset")

    top_ha <- ComplexHeatmap::HeatmapAnnotation(
        df = as.data.frame(annotation_list),
        col = color_list,
        show_legend = TRUE,
        show_annotation_name = FALSE,
        annotation_legend_param = legend_params
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
            cluster_columns = if(pseudo_bulk && cluster_cols) TRUE else FALSE,
            column_order = if(pseudo_bulk && cluster_cols) NULL else cell_subset[["cell_id"]],
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
            row_names_gp = row_gp,
            cluster_rows = TRUE,
            show_column_names = FALSE,
            cluster_columns = if(pseudo_bulk && cluster_cols) TRUE else FALSE,
            column_order = if(pseudo_bulk && cluster_cols) NULL else cell_subset[["cell_id"]],
            top_annotation = top_ha,
            border = TRUE,
            use_raster = TRUE,
            raster_quality = 2
        )
    }

    return(ht)
}

#' @title Plot Barplots for Top Loading Gene Fold Changes (Pseudo-Bulk)
#'
#' @description
#' This internal helper function creates a barplot visualization showing pseudo-bulk fold changes
#' between different cell groups for top loading genes.
#'
#' @details
#' This function generates a ggplot2 barplot where genes are arranged vertically
#' in hierarchically clustered order (identical to heatmap gene ordering), and horizontal bars show log2 fold changes
#' for different pseudo-bulk comparisons vs reference. The function creates an internal heatmap
#' object with the same parameters to extract the exact gene clustering order, ensuring
#' consistency between plot types.
#'
#' Bar colors:
#' \itemize{
#'   \item Green: Query non-anomaly (pseudo-bulk) vs Reference (pseudo-bulk)
#'   \item Yellow: All Query (pseudo-bulk) vs Reference (pseudo-bulk) (shown when show_all_query = TRUE and anomaly data available)
#'   \item Red: Query anomaly (pseudo-bulk) vs Reference (pseudo-bulk)
#' }
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject} containing
#'          expression data and analysis results.
#' @param cell_type A character string specifying the cell type to visualize.
#' @param available_pcs A character vector of principal components to include in analysis.
#' @param plot_by A character string indicating gene selection criterion ("top_loading" or "p_adjusted").
#' @param n_genes An integer specifying the number of top genes to display per PC.
#' @param significance_threshold A numeric value between 0 and 1 for significance annotation.
#' @param show_anomalies Logical indicating whether to show anomaly-related bars.
#' @param show_all_query Logical indicating whether to show the yellow bar for all query vs reference
#'                       comparison. Default is TRUE. When FALSE, only green and red bars are shown.
#'
#' @keywords internal
#'
#' @return A ggplot2 object ready for display, or NULL if no genes meet selection criteria.
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
plotBarplot <- function(x, cell_type, available_pcs, plot_by,
                        n_genes, significance_threshold, show_anomalies,
                        show_all_query = TRUE) {

    # Validate required parameters
    if (is.null(n_genes)) {
        stop("'plot_type = \"barplot\"' requires 'n_genes' to be set.")
    }

    # Validate show_all_query parameter
    if (!is.logical(show_all_query) || length(show_all_query) != 1) {
        stop("show_all_query must be a logical value.")
    }

    # Create internal heatmap object to extract gene order
    # Use the same parameters as would be used for heatmap plotting
    internal_heatmap <- tryCatch({
        plotHeatmap(x, cell_type, available_pcs, plot_by, n_genes,
                    significance_threshold, show_anomalies,
                    pseudo_bulk = TRUE, cluster_cols = FALSE)
    }, error = function(e) {
        stop("Failed to create internal heatmap for gene ordering: ", e$message)
    })

    if (is.null(internal_heatmap)) {
        warning("No genes met the selection criteria for the barplot.")
        return(NULL)
    }

    # Extract gene order from the heatmap
    gene_order_clustered <- extractGeneOrder(internal_heatmap)

    # Get expression data for these genes
    cell_subset <- x[["cell_metadata"]][x[["cell_metadata"]][["cell_type"]] == cell_type, ]
    expr_matrix <- as.matrix(
        x[["expression_data"]][gene_order_clustered, cell_subset[["cell_id"]], drop = FALSE]
    )

    # Check if anomaly data is available
    has_anomaly_data <- show_anomalies && "anomaly_status" %in% names(cell_subset)

    # Create pseudo-bulk profiles for each group
    if (has_anomaly_data) {
        # Only show anomalies for query cells, treat all reference as normal
        cell_subset_modified <- cell_subset
        cell_subset_modified[["anomaly_status"]][cell_subset_modified[["dataset"]] == "Reference"] <- "Normal"

        # 3 categories: Query_Normal, Query_Anomaly, Reference_Normal
        cell_subset_modified[["group"]] <- paste(cell_subset_modified[["dataset"]],
                                                 cell_subset_modified[["anomaly_status"]],
                                                 sep = "_")
        group_order <- c("Query_Normal", "Query_Anomaly", "Reference_Normal")
        cell_subset <- cell_subset_modified
    } else {
        # 2 categories: Reference, Query
        cell_subset[["group"]] <- cell_subset[["dataset"]]
        group_order <- c("Reference", "Query")
    }

    # Calculate pseudo-bulk expression profiles
    pseudo_bulk_matrix <- matrix(0, nrow = nrow(expr_matrix), ncol = length(group_order))
    rownames(pseudo_bulk_matrix) <- rownames(expr_matrix)
    colnames(pseudo_bulk_matrix) <- group_order

    for (group in group_order) {
        group_cells <- cell_subset[cell_subset[["group"]] == group, "cell_id"]
        if (length(group_cells) > 0) {
            if (length(group_cells) == 1) {
                pseudo_bulk_matrix[, group] <- expr_matrix[, group_cells]
            } else {
                pseudo_bulk_matrix[, group] <- rowMeans(expr_matrix[, group_cells, drop = FALSE])
            }
        }
    }

    # Calculate fold changes for different comparisons
    fold_change_data <- list()

    # Reference pseudo-bulk (always Reference_Normal or Reference)
    ref_column <- if (has_anomaly_data) "Reference_Normal" else "Reference"

    for (gene in gene_order_clustered) {
        ref_expr <- pseudo_bulk_matrix[gene, ref_column]

        # Query non-anomaly vs Reference (always calculated)
        query_normal_column <- if (has_anomaly_data) "Query_Normal" else "Query"
        query_normal_expr <- pseudo_bulk_matrix[gene, query_normal_column]
        query_normal_fc <- query_normal_expr - ref_expr  # Log2 fold change

        # All Query vs Reference (calculated when anomaly data available AND show_all_query is TRUE)
        if (has_anomaly_data && show_all_query) {
            # Calculate overall query mean from both normal and anomaly
            query_normal_cells <- cell_subset[cell_subset[["group"]] == "Query_Normal", "cell_id"]
            query_anomaly_cells <- cell_subset[cell_subset[["group"]] == "Query_Anomaly", "cell_id"]
            all_query_cells <- c(query_normal_cells, query_anomaly_cells)

            if (length(all_query_cells) > 0) {
                query_all_expr <- mean(expr_matrix[gene, all_query_cells], na.rm = TRUE)
                query_all_fc <- query_all_expr - ref_expr
            } else {
                query_all_fc <- NA
            }
        } else {
            query_all_fc <- NA
        }

        # Query anomaly vs Reference (calculated when anomaly data available)
        if (has_anomaly_data) {
            query_anomaly_expr <- pseudo_bulk_matrix[gene, "Query_Anomaly"]
            query_anomaly_fc <- query_anomaly_expr - ref_expr
        } else {
            query_anomaly_fc <- NA
        }

        # Store fold change data
        gene_fc_data <- data.frame(
            gene = gene,
            query_normal_fc = query_normal_fc,
            query_all_fc = query_all_fc,
            query_anomaly_fc = query_anomaly_fc,
            stringsAsFactors = FALSE
        )

        fold_change_data[[length(fold_change_data) + 1]] <- gene_fc_data
    }

    # Combine fold change data
    fc_df <- do.call(rbind, fold_change_data)

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

    # Reshape data for plotting
    plot_data_list <- list()

    # Query non-anomaly vs Reference (always shown)
    normal_data <- data.frame(
        gene = fc_df[["gene"]],
        fold_change = fc_df[["query_normal_fc"]],
        comparison = "Query Non-Anomaly vs Reference",
        color_group = "normal",
        stringsAsFactors = FALSE
    )
    plot_data_list[[length(plot_data_list) + 1]] <- normal_data

    # All Query vs Reference (shown when anomaly data available AND show_all_query is TRUE)
    if (has_anomaly_data && show_all_query && !all(is.na(fc_df[["query_all_fc"]]))) {
        all_data <- data.frame(
            gene = fc_df[["gene"]],
            fold_change = fc_df[["query_all_fc"]],
            comparison = "All Query vs Reference",
            color_group = "all",
            stringsAsFactors = FALSE
        )
        plot_data_list[[length(plot_data_list) + 1]] <- all_data
    }

    # Query anomaly vs Reference (shown when anomaly data available)
    if (has_anomaly_data && !all(is.na(fc_df[["query_anomaly_fc"]]))) {
        anomaly_data <- data.frame(
            gene = fc_df[["gene"]],
            fold_change = fc_df[["query_anomaly_fc"]],
            comparison = "Query Anomaly vs Reference",
            color_group = "anomaly",
            stringsAsFactors = FALSE
        )
        plot_data_list[[length(plot_data_list) + 1]] <- anomaly_data
    }

    # Combine plot data
    if (length(plot_data_list) == 0) {
        warning("No fold change data available for plotting.")
        return(NULL)
    }

    plot_df <- do.call(rbind, plot_data_list)
    plot_df <- plot_df[!is.na(plot_df[["fold_change"]]), ]

    if (nrow(plot_df) == 0) {
        warning("No valid fold change data for plotting.")
        return(NULL)
    }

    # Set gene factor levels in clustered order (reversed for y-axis display)
    plot_df[["gene"]] <- factor(plot_df[["gene"]], levels = rev(gene_order_clustered))

    # Set comparison factor levels for proper ordering
    comparison_levels <- c("All Query vs Reference",
                           "Query Non-Anomaly vs Reference",
                           "Query Anomaly vs Reference")
    plot_df[["comparison"]] <- factor(plot_df[["comparison"]], levels = comparison_levels)

    # Set factor levels on color_group to control bar order
    plot_df[["color_group"]] <- factor(plot_df[["color_group"]],
                                       levels = c("anomaly", "normal", "all"))

    bar_colors <- c("normal" = "#9E9E9E",
                    "all" = "#666666",
                    "anomaly" = "#DC2F41")

    # Prepare gene labels with significance indicators
    gene_labels <- sapply(gene_order_clustered, function(g) {
        paste0(g, ifelse(g %in% significant_genes, "*", ""))
    })
    names(gene_labels) <- gene_order_clustered

    # Create the plot
    p <- ggplot2::ggplot(data = plot_df,
                         ggplot2::aes(x = .data[["fold_change"]],
                                      y = .data[["gene"]],
                                      fill = .data[["color_group"]])) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                          alpha = 0.8, width = 0.4) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            color = "grey50", alpha = 0.7) +
        ggplot2::scale_fill_manual(name = "Comparison",
                                   values = bar_colors,
                                   labels = c("normal" = "Query Non-Anomaly vs Ref",
                                              "all" = "All Query vs Ref",
                                              "anomaly" = "Query Anomaly vs Ref"),
                                   guide = ggplot2::guide_legend(reverse = TRUE)) +
        ggplot2::scale_y_discrete(labels = rev(gene_labels), name = NULL) +
        ggplot2::scale_x_continuous(name = "Log2 Fold Change (Pseudo-Bulk)",
                                    expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
        ggplot2::labs(title = NULL, subtitle = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(face = "bold", size = 10),
            axis.text.x = ggplot2::element_text(size = 9),
            axis.title.x = ggplot2::element_text(face = "bold", size = 11),
            legend.position = "right",
            legend.title = ggplot2::element_text(face = "bold", size = 10),
            legend.text = ggplot2::element_text(size = 9),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
            panel.border = ggplot2::element_rect(color = "grey70", linewidth = 0.5)
        )

    return(p)
}

#' @title Plot Boxplots for Top Loading Gene Shifts (Two-Panel Summary using ggplot2)
#'
#' @description
#' This internal helper function creates a comprehensive two-panel summary plot
#' displaying gene expression distributions and principal component loadings.
#' The visualization uses ggplot2 faceting to create side-by-side panels.
#' Optionally includes anomaly status information using visual cues.
#'
#' @details
#' This function generates a dual-panel ggplot2 visualization where the left panel
#' shows horizontal boxplots of gene expression distributions comparing Reference
#' and Query datasets, while the right panel displays PC loading values as points
#' with adjusted p-values. Gene selection is based on the union of top genes
#' across specified principal components.
#'
#' When anomaly information is available and requested, anomalous cells are
#' distinguished using dashed boxplot borders (normal cells have solid borders).
#' This approach avoids color conflicts between dataset identification (fill colors)
#' and PC identification (point colors/shapes).
#'
#' Visual encoding:
#' \itemize{
#'   \item Dataset: Fill colors (Reference = blue, Query = red)
#'   \item Anomaly status: Line types (Normal = solid, Anomaly = dashed borders)
#'   \item PC identity: Point colors and shapes in loading panel
#' }
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject} containing
#'   expression data and analysis results.
#' @param cell_type A character string specifying the cell type to visualize.
#' @param available_pcs A character vector of principal components to include in analysis.
#' @param plot_by A character string indicating gene selection criterion ("top_loading" or "p_adjusted").
#' @param n_genes An integer specifying the number of top genes to display per PC.
#' @param significance_threshold A numeric value between 0 and 1 for significance annotation.
#' @param show_anomalies Logical indicating whether to show anomaly annotations using
#'   line type differences in boxplot borders.
#'
#' @keywords internal
#'
#' @return A ggplot2 object ready for display, or NULL if no genes meet selection criteria.
#'   The returned plot contains:
#'   \itemize{
#'     \item Left panel: Expression boxplots with dataset-specific fill colors
#'     \item Right panel: PC loading scatter points with PC-specific colors/shapes
#'     \item Gene labels on y-axis with significance indicators (*)
#'     \item P-value annotations on secondary y-axis
#'     \item Legend showing dataset, PC information, and anomaly status (if applicable)
#'   }
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plotHeatmap}}, \code{\link{plot.calculateTopLoadingGeneShiftsObject}}
#'
# Helper boxplot function
plotBoxplot <- function(x, cell_type, available_pcs, plot_by,
                        n_genes, significance_threshold, show_anomalies) {

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

    # Create expression data frame with cell-level metadata
    expr_plot_list <- lapply(final_gene_data[["gene"]], function(g) {
        gene_data <- data.frame(
            gene = g,
            value = as.numeric(expr_subset[g, ]),
            dataset = cell_subset_meta[["dataset"]],
            metric = "Expression",
            pc_group = NA_character_,
            stringsAsFactors = FALSE
        )

        # Add anomaly status if available and requested
        if (show_anomalies && "anomaly_status" %in% names(cell_subset_meta)) {
            gene_data[["anomaly_status"]] <- cell_subset_meta[["anomaly_status"]]
        }

        return(gene_data)
    })
    expr_plot_df <- do.call(rbind, expr_plot_list)

    # Prepare loading data for plotting (PC-specific information)
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

    # Add anomaly status column to loading data (not applicable for PC loadings)
    if (show_anomalies && "anomaly_status" %in% names(cell_subset_meta)) {
        loading_plot_df[["anomaly_status"]] <- NA_character_
    }

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

    # Format p-value labels for secondary axis
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
    pc_shape_palette <- c(16, 17, 15, 18, 8)  # circle, triangle, square, diamond, star

    pc_colors <- pc_color_palette[1:num_pcs]
    pc_shapes <- pc_shape_palette[1:num_pcs]
    names(pc_colors) <- available_pcs
    names(pc_shapes) <- available_pcs

    # Create facet labeller for panel titles
    facet_labeller <- ggplot2::labeller(
        metric = c("Expression" = "Log-Normalized Expression",
                   "PC Loading" = "PC Loading")
    )

    # Build the base ggplot object
    p <- ggplot2::ggplot(
        data = full_plot_df,
        mapping = ggplot2::aes(y = .data[["gene"]], x = .data[["value"]])
    )

    # Add expression panel elements with conditional anomaly visualization
    if (show_anomalies && "anomaly_status" %in% names(full_plot_df)) {
        # Boxplots with anomaly-aware border styling using linetype
        p <- p + ggplot2::geom_boxplot(
            data = ~subset(., metric == "Expression"),
            ggplot2::aes(fill = .data[["dataset"]],
                         linetype = .data[["anomaly_status"]]),
            alpha = 0.8,
            outlier.size = 0.5,
            size = 0.8  # Slightly thicker borders for better linetype visibility
        ) +
            ggplot2::scale_linetype_manual(
                name = "Status",
                values = c("Normal" = "solid", "Anomaly" = "dashed"),
                na.value = "solid",
                guide = ggplot2::guide_legend(
                    order = 3,
                    override.aes = list(color = "black", fill = "grey80")
                )
            )
    } else {
        # Standard boxplots without anomaly information
        p <- p + ggplot2::geom_boxplot(
            data = ~subset(., metric == "Expression"),
            ggplot2::aes(fill = .data[["dataset"]]),
            alpha = 0.8,
            outlier.size = 0.5
        )
    }

    # Add PC loading panel elements
    p <- p +
        # Reference line at zero for PC loadings
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed",
            color = "grey50",
            alpha = 0.7
        ) +
        # PC loading points with color and shape encoding
        ggplot2::geom_point(
            data = ~subset(., metric == "PC Loading"),
            ggplot2::aes(color = .data[["pc_group"]],
                         shape = .data[["pc_group"]]),
            size = 3,
            alpha = 0.8
        ) +

        # Panel configuration with custom labelling
        ggplot2::facet_grid(
            ~ metric,
            scales = "free_x",
            switch = "x",
            labeller = facet_labeller
        ) +

        # Define all aesthetic scales
        ggplot2::scale_fill_manual(
            name = "Dataset",
            values = dataset_colors,
            guide = ggplot2::guide_legend(
                order = 1,
                override.aes = list(linetype = "solid")
            )
        ) +
        ggplot2::scale_color_manual(
            name = "PC",
            values = pc_colors,
            labels = pc_legend_labels,
            guide = ggplot2::guide_legend(order = 2)
        ) +
        ggplot2::scale_shape_manual(
            name = "PC",
            values = pc_shapes,
            labels = pc_legend_labels,
            guide = ggplot2::guide_legend(order = 2)
        ) +
        ggplot2::scale_y_discrete(
            labels = rev(gene_labels),
            name = NULL
        ) +
        ggplot2::scale_x_continuous(
            sec.axis = ggplot2::dup_axis(
                name = "Adj. P-value",
                breaks = NULL,
                labels = rev(p_value_labels)
            ),
            name = NULL,
            expand = ggplot2::expansion(mult = c(0.05, 0.15))
        ) +

        # Labels and legend configuration
        ggplot2::labs(
            title = NULL,
            subtitle = NULL
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(order = 1),
            color = ggplot2::guide_legend(order = 2),
            shape = ggplot2::guide_legend(order = 2),
            linetype = if(show_anomalies && "anomaly_status" %in% names(full_plot_df)) {
                ggplot2::guide_legend(order = 3)
            } else {
                "none"
            }
        ) +

        # Theme customization for publication-ready appearance
        ggplot2::theme_bw() +
        ggplot2::theme(
            # Remove x-axis title (shown via facet labels)
            axis.title.x = ggplot2::element_blank(),
            # Bold gene names for emphasis
            axis.text.y.left = ggplot2::element_text(face = "bold", size = 10),
            # P-value labels styling
            axis.text.y.right = ggplot2::element_text(size = 8),
            # Legend positioning and styling
            legend.position = "right",
            legend.box = "vertical",
            legend.key.size = ggplot2::unit(0.8, "lines"),
            legend.text = ggplot2::element_text(size = 9),
            legend.title = ggplot2::element_text(size = 10, face = "bold"),
            # Panel styling
            strip.background = ggplot2::element_blank(),
            strip.placement = "outside",
            strip.text.x = ggplot2::element_text(face = "bold", size = 11),
            # Grid and panel borders
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "grey70", linewidth = 0.5)
        )

    return(p)
}

#' @title Extract Gene Order from ComplexHeatmap Object
#'
#' @description
#' Extracts the hierarchically clustered gene order from a ComplexHeatmap Heatmap object.
#'
#' @details
#' This function initializes a ComplexHeatmap object (without displaying it) to perform
#' hierarchical clustering, then extracts the resulting gene order. This allows matching
#' gene ordering between heatmap and other plot types. The function uses a null graphics
#' device to initialize the heatmap clustering without actually drawing the plot.
#'
#' @param heatmap_object A ComplexHeatmap Heatmap object (undrawn or drawn).
#'
#' @keywords internal
#'
#' @return A character vector of gene names in the clustered order (top to bottom as they appear in heatmap).
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @importFrom grDevices dev.off pdf
#'
# Function to extract genes used by heatmap plot
extractGeneOrder <- function(heatmap_object) {

    # Check if it's a ComplexHeatmap object
    if (!inherits(heatmap_object, "Heatmap")) {
        stop("Input must be a ComplexHeatmap Heatmap object")
    }

    tryCatch({
        # Open a null device to draw without displaying
        pdf(NULL)
        on.exit(dev.off(), add = TRUE)

        # Draw to initialize the heatmap clustering
        drawn_ht <- ComplexHeatmap::draw(heatmap_object, show_heatmap_legend = FALSE,
                                         show_annotation_legend = FALSE)

        # Extract row order from the drawn heatmap
        gene_order_indices <- ComplexHeatmap::row_order(drawn_ht)

        # Get the matrix and extract gene names in clustered order
        heatmap_matrix <- heatmap_object@matrix
        if (is.null(rownames(heatmap_matrix))) {
            stop("Heatmap matrix has no row names (gene names)")
        }

        # Return genes in clustered order
        clustered_genes <- rownames(heatmap_matrix)[gene_order_indices]
        return(clustered_genes)

    }, error = function(e) {
        stop("Failed to extract gene order from heatmap object. Error: ", e$message)
    })
}
