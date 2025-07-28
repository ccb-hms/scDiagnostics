#' @title Plot Top Loading Gene Expression Shifts
#'
#' @description
#' This function creates boxplots showing expression distributions for top loading genes
#' that exhibit distributional differences between query and reference datasets.
#'
#' @details
#' This function visualizes the results from \code{calculateTopLoadingGeneShifts} by creating
#' faceted boxplots. Each facet represents a principal component with its variance explained.
#' Within each facet, boxplots show expression distributions for selected genes, comparing
#' query and reference datasets for a specific cell type. Genes can be selected either by
#' highest absolute loadings or most significant p-adjusted values.
#'
#' @param x An object of class \code{calculateTopLoadingGeneShiftsObject} containing the output of the \code{calculateTopLoadingGeneShifts} function.
#' @param cell_type A character string specifying the cell type to plot (must be exactly one).
#' @param pc_subset A numeric vector specifying which principal components to plot. Default is 1:3.
#' @param plot_by A character string specifying gene selection method. Either "top_loading" or "p_adjusted". Default is "p_adjusted".
#' @param n_genes Number of genes to show per PC. Default is 10.
#' @param significance_threshold P-adjusted threshold for significance annotation. Default is 0.05.
#' @param ... Additional arguments to be passed to the plotting functions.
#'
#' @return A ggplot object showing boxplots of gene expression by PC and dataset.
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
# Function to plot results of calculateTopLoadingGeneShifts
plot.calculateTopLoadingGeneShiftsObject <- function(x,
                                                     cell_type,
                                                     pc_subset = 1:3,
                                                     plot_by = c("p_adjusted", "top_loading"),
                                                     n_genes = 10,
                                                     significance_threshold = 0.05,
                                                     ...) {

    # Input validation
    if (missing(cell_type) || length(cell_type) != 1) {
        stop("cell_type must be specified and be exactly one cell type")
    }

    if (!is.numeric(n_genes) || n_genes <= 0) {
        stop("n_genes must be a positive integer")
    }

    # Match plot_by argument
    plot_by <- match.arg(plot_by)

    # Check if x has the expected structure
    required_elements <- c("expression_data", "cell_metadata", "gene_metadata", "percent_var")
    if (!all(required_elements %in% names(x))) {
        stop("x must contain expression_data, cell_metadata, gene_metadata, and percent_var elements")
    }

    # Get PC names in proper order (PC1, PC2, ... not PC1, PC10, PC2)
    pc_names <- paste0("PC", pc_subset)
    available_pcs <- intersect(pc_names, names(x))

    if (length(available_pcs) == 0) {
        stop("No requested PCs found in x")
    }

    # Reorder to match requested order
    available_pcs <- pc_names[pc_names %in% available_pcs]

    # Check if cell type exists
    if (!cell_type %in% x[["cell_metadata"]][["cell_type"]]) {
        stop(paste("Cell type", cell_type, "not found in results"))
    }

    # Create PC labels with variance explained from stored percent_var
    pc_labels <- paste0("PC", pc_subset, " (", sprintf("%.1f%%", x[["percent_var"]][paste0("PC", pc_subset)]), ")")
    names(pc_labels) <- paste0("PC", pc_subset)

    # Filter to available PCs
    pc_labels <- pc_labels[available_pcs]

    # Collect data for plotting
    plot_data_list <- list()
    gene_annotation_list <- list()

    for (pc_name in available_pcs) {
        pc_results <- x[[pc_name]]

        if (nrow(pc_results) == 0) {
            next
        }

        # Filter for the specified cell type
        pc_cell_data <- pc_results[pc_results[["cell_type"]] == cell_type, ]

        if (nrow(pc_cell_data) == 0) {
            warning(paste("No data found for cell type", cell_type, "in", pc_name))
            next
        }

        # Select genes based on plot_by criterion
        if (plot_by == "top_loading") {
            # Order by absolute loading (descending)
            gene_order <- order(abs(pc_cell_data[["loading"]]), decreasing = TRUE)
        } else if (plot_by == "p_adjusted") {
            # Order by p_adjusted (ascending - most significant first)
            gene_order <- order(pc_cell_data[["p_adjusted"]])
        }

        # Get top n_genes
        selected_indices <- gene_order[1:min(n_genes, length(gene_order))]
        selected_data <- pc_cell_data[selected_indices, ]

        # Get expression data for selected genes
        selected_genes <- selected_data[["gene"]]

        # Filter cell metadata for this cell type
        cell_subset <- x[["cell_metadata"]][["cell_type"]] == cell_type
        relevant_cells <- x[["cell_metadata"]][cell_subset, ]

        # Extract expression data
        expr_subset <- x[["expression_data"]][selected_genes,
                                              relevant_cells[["cell_id"]],
                                              drop = FALSE]

        # Create long format data for plotting
        for (gene in selected_genes) {
            gene_expr <- expr_subset[gene, ]

            gene_plot_data <- data.frame(
                PC = pc_labels[pc_name],  # Use labeled PC name
                Gene = gene,
                Expression = as.numeric(gene_expr),
                Dataset = relevant_cells[["dataset"]],
                row.names = NULL,  # Explicitly set row.names to NULL to avoid warnings
                stringsAsFactors = FALSE
            )

            plot_data_list[[length(plot_data_list) + 1]] <- gene_plot_data

            # Store gene annotation info
            gene_info <- selected_data[selected_data[["gene"]] == gene, ]
            gene_annotation_list[[length(gene_annotation_list) + 1]] <- data.frame(
                PC = pc_labels[pc_name],  # Use labeled PC name
                Gene = gene,
                p_adjusted = gene_info[["p_adjusted"]],
                loading = gene_info[["loading"]],
                significant = gene_info[["significant"]],
                row.names = NULL,  # Explicitly set row.names to NULL to avoid warnings
                stringsAsFactors = FALSE
            )
        }
    }

    # Combine all data
    if (length(plot_data_list) == 0) {
        stop("No data available for plotting with the specified parameters")
    }

    plot_data <- do.call(rbind, plot_data_list)
    gene_annotations <- do.call(rbind, gene_annotation_list)

    # Set proper factor levels for PC ordering (with labels)
    pc_label_levels <- pc_labels[available_pcs]
    plot_data[["PC"]] <- factor(plot_data[["PC"]], levels = pc_label_levels)
    gene_annotations[["PC"]] <- factor(gene_annotations[["PC"]], levels = pc_label_levels)

    # Set Dataset factor levels to control boxplot order (Reference left, Query right)
    plot_data[["Dataset"]] <- factor(plot_data[["Dataset"]], levels = c("Reference", "Query"))

    # Create a Gene_PC interaction variable for proper ordering within each facet
    # This ensures each PC has its own gene ordering
    if (plot_by == "top_loading") {
        # For top_loading: order by absolute loading (descending) within each PC
        gene_annotations[["sort_value"]] <- -abs(gene_annotations[["loading"]])
    } else {
        # For p_adjusted: order by p_adjusted (ascending) within each PC
        gene_annotations[["sort_value"]] <- gene_annotations[["p_adjusted"]]
    }

    # Create a unique identifier for each gene-PC combination for proper ordering
    gene_annotations[["Gene_PC"]] <- paste(gene_annotations[["Gene"]],
                                           gene_annotations[["PC"]],
                                           sep = "_")

    # Order within each PC
    gene_annotations <- gene_annotations[order(gene_annotations[["PC"]],
                                               gene_annotations[["sort_value"]]), ]

    # Create the same Gene_PC variable for plot_data
    plot_data[["Gene_PC"]] <- paste(plot_data[["Gene"]],
                                    plot_data[["PC"]],
                                    sep = "_")

    # Set factor levels for Gene_PC to maintain proper ordering
    plot_data[["Gene_PC"]] <- factor(plot_data[["Gene_PC"]],
                                     levels = unique(gene_annotations[["Gene_PC"]]))
    gene_annotations[["Gene_PC"]] <- factor(gene_annotations[["Gene_PC"]],
                                            levels = unique(gene_annotations[["Gene_PC"]]))

    # Create better significance labels
    gene_annotations[["p_formatted"]] <- ifelse(
        gene_annotations[["p_adjusted"]] < 0.001,
        "p < 0.001",  # Changed format for very small p-values
        paste0("p = ", sprintf("%.3f", gene_annotations[["p_adjusted"]]))  # Added "p = " prefix
    )

    # Calculate number of columns for facet wrap (max 3 per row)
    n_facets <- length(available_pcs)
    n_cols <- min(n_facets, 3)

    # Create the main plot using Gene_PC for proper ordering
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["Gene_PC"]],
                                                 y = .data[["Expression"]],
                                                 fill = .data[["Dataset"]])) +
        ggplot2::geom_boxplot(alpha = 0.7,
                              outlier.size = 0.5,  # Show outliers with small size
                              width = 0.7) +
        ggplot2::facet_wrap(~ PC, scales = "free",  # Free both x and y scales
                            ncol = n_cols) +  # Dynamic number of columns
        ggplot2::scale_fill_manual(values = c("Reference" = "#3498DB", "Query" = "#E74C3C"),
                                   name = "Dataset") +
        ggplot2::scale_x_discrete(labels = function(x) {
            # Extract just the gene name from Gene_PC for x-axis labels
            sapply(strsplit(x, "_"), function(parts) parts[1])
        }) +
        ggplot2::labs(
            title = paste0("Gene Expression Shifts for ", cell_type),
            subtitle = paste0("Genes selected by ",
                              ifelse(plot_by == "top_loading", "highest absolute loadings",
                                     "most significant adjusted p-values")),
            x = "",  # No x-axis label
            y = ""   # No y-axis label
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = "white",
                                                     color = "black",
                                                     linewidth = 0.5),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold",
                                               hjust = 0),  # Left aligned title
            plot.subtitle = ggplot2::element_text(size = 11,
                                                  hjust = 0,  # Left aligned subtitle
                                                  color = "gray30"),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(angle = 45,
                                                hjust = 1,
                                                size = 9),
            legend.position = "right",  # Legend on right side
            legend.title = ggplot2::element_text(size = 11),
            legend.text = ggplot2::element_text(size = 10),
            strip.text = ggplot2::element_text(size = 11, face = "plain")  # No bold facet titles
        )

    # Add significance annotations with better positioning
    # Calculate y position for annotations based on data range per facet
    y_ranges <- aggregate(Expression ~ PC, data = plot_data,
                          FUN = function(x) c(min = min(x, na.rm = TRUE),
                                              max = max(x, na.rm = TRUE)))

    # Create a data frame for y positions
    y_positions <- data.frame(
        PC = y_ranges[["PC"]],
        y_min = y_ranges[["Expression"]][, "min"],
        y_max = y_ranges[["Expression"]][, "max"]
    )

    # Merge with gene annotations
    gene_annotations <- merge(gene_annotations, y_positions, by = "PC")
    gene_annotations[["y_annotation"]] <- gene_annotations[["y_max"]] +
        0.05 * (gene_annotations[["y_max"]] - gene_annotations[["y_min"]])

    # Add p-values (no asterisks)
    p <- p + ggplot2::geom_text(
        data = gene_annotations,
        ggplot2::aes(x = .data[["Gene_PC"]], y = .data[["y_annotation"]],
                     label = .data[["p_formatted"]]),
        inherit.aes = FALSE,
        size = 2.8,
        color = ifelse(gene_annotations[["significant"]], "#D32F2F", "gray50"),
        vjust = 0,
        angle = 0
    )

    return(p)
}
