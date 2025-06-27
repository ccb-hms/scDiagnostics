#' @title Visualize gene expression on a dimensional reduction plot
#'
#' @description
#' This function plots gene expression on a dimensional reduction plot using methods like t-SNE, UMAP, or PCA. Each single cell is color-coded based on the expression of a specific gene or feature.
#'
#' @param se_object An object of class \code{\linkS4class{SingleCellExperiment}} containing log-transformed expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param method The reduction method to use for visualization. It should be one of the supported methods: "TSNE", "UMAP", or "PCA".
#' @param pc_subset An optional vector specifying the principal components (PCs) to include in the plot if method = "PCA".
#'        Default is 1:5.
#' @param feature A character string representing the name of the gene or feature to be visualized.
#' @param cell_type_col The column name in the \code{colData} of \code{se_object} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#'
#' @importFrom SummarizedExperiment assay
#' @import SingleCellExperiment
#'
#' @return A ggplot object representing the dimensional reduction plot with gene expression.
#'
#' @export
#'
#' @examples
#' # Load data
#' data("query_data")
#'
#' # Plot gene expression on PCA plot
#' plotGeneExpressionDimred(se_object = query_data,
#'                          cell_type_col = "SingleR_annotation",
#'                          method = "PCA",
#'                          pc_subset = 1:5,
#'                          feature = "CD8A",
#'                          cell_types = "CD4")
#'
# Function to plot the gene expression of a gene for one or more cell types
plotGeneExpressionDimred <- function(se_object,
                                     method = c("TSNE", "UMAP", "PCA"),
                                     pc_subset = 1:5,
                                     feature,
                                     cell_type_col = NULL,
                                     cell_types = NULL,
                                     assay_name = "logcounts") {

    # Check standard input arguments
    argumentCheck(query_data = se_object,
                  pc_subset_query = pc_subset,
                  assay_name = assay_name,
                  query_cell_type_col = cell_type_col,
                  cell_types = cell_types)

    # Match arguments
    method <- match.arg(method)

    # Check if feature is available
    if (!feature %in% rownames(assay(se_object, assay_name))) {
        stop("Specified feature does not exist in the expression matrix.")
    }

    # Handle cell type filtering
    filtered_to_specific_types <- FALSE
    if (!is.null(cell_type_col)) {
        # Check if cell type column exists
        if (!cell_type_col %in% colnames(colData(se_object))) {
            stop("Specified cell_type_col does not exist in colData.")
        }

        # Get available cell types
        available_cell_types <- unique(colData(se_object)[[cell_type_col]])

        # If specific cell types are requested, filter to those
        if (!is.null(cell_types)) {
            # Check if requested cell types exist
            missing_types <- setdiff(cell_types, available_cell_types)
            if (length(missing_types) > 0) {
                warning("The following cell types were not found: ",
                        paste(missing_types, collapse = ", "))
            }

            # Filter to cells of specified types
            cell_mask <- colData(se_object)[[cell_type_col]] %in% cell_types
            filtered_to_specific_types <- TRUE
        } else {
            # Use all cells if no specific types requested
            cell_mask <- rep(TRUE, ncol(se_object))
            cell_types <- available_cell_types
        }

        # Filter the object
        se_object <- se_object[, cell_mask]
    }

    # Extract gene expression vector (after potential filtering)
    expression <- assay(se_object, assay_name)[feature, ]

    if(method %in% c("TSNE", "UMAP")){

        # Extract dimension reduction coordinates from SingleCellExperiment object
        reduction <- reducedDim(se_object, method)

        # Prepare data for plotting
        df <- data.frame(Dim1 = reduction[, 1], Dim2 = reduction[, 2],
                         Expression = expression)

        # Add cell type information if available
        if (!is.null(cell_type_col)) {
            df$CellType <- colData(se_object)[[cell_type_col]]
        }

        # Create the plot object with better color gradient
        plot_obj <- ggplot2::ggplot(df, ggplot2::aes(
            x = .data[["Dim1"]],
            y = .data[["Dim2"]])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["Expression"]]),
                                alpha = 0.7, size = 1.2) +
            ggplot2::scale_color_gradient(low = "lightgray", high = "red",
                                          name = paste(feature, "\nExpression")) +
            ggplot2::xlab("Dimension 1") +
            ggplot2::ylab("Dimension 2") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                strip.background = ggplot2::element_rect(fill = "grey85",
                                                         color = "grey70"),
                strip.text = ggplot2::element_text(size = 10,
                                                   face = "bold",
                                                   color = "black"),
                axis.title = ggplot2::element_text(size = 12),
                axis.text = ggplot2::element_text(size = 10),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white",
                                                         color = "black"),
                legend.position = "right",
                plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
                plot.background = ggplot2::element_rect(fill = "white"))

        # Add title with cell type info if filtered
        if (!is.null(cell_type_col) && !is.null(cell_types)) {
            plot_obj <- plot_obj +
                ggplot2::ggtitle(paste0(feature, " Expression in ",
                                        paste(cell_types, collapse = ", "),
                                        " Cells"))
        }

    } else if (method == "PCA"){

        # Check input for pc_subset
        if(!all(pc_subset %in% seq_len(ncol(reducedDim(se_object, "PCA")))))
            stop("\'pc_subset\' is out of range.")

        # PCA data
        plot_mat <- reducedDim(se_object, "PCA")[, pc_subset]

        # Create PC column names - include variance explained only if using all cell types
        if (filtered_to_specific_types) {
            plot_names <- paste0("PC", pc_subset)
        } else {
            plot_names <- paste0(
                "PC", pc_subset, " (",
                sprintf("%.1f%%", attributes(
                    reducedDim(se_object, "PCA"))[["percentVar"]][pc_subset]), ")")
        }

        # Create a new data frame with selected PCs
        pc_df <- data.frame(matrix(0, nrow = nrow(plot_mat),
                                   ncol = length(pc_subset)))
        colnames(pc_df) <- plot_names

        for (i in 1:length(pc_subset)) {
            pc_df[, i] <- plot_mat[, i]
        }

        # Add expression data
        pc_df[["Expression"]] <- expression

        # Add cell type information if available
        if (!is.null(cell_type_col)) {
            pc_df[["CellType"]] <- colData(se_object)[[cell_type_col]]
        }

        # Create a simple plot to extract the legend with better color gradient
        legend_plot <- ggplot2::ggplot(pc_df,
                                       ggplot2::aes(x = pc_df[,1],
                                                    y = pc_df[,2])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["Expression"]])) +
            ggplot2::scale_color_gradient(low = "gray", high = "blue",
                                          name = paste(feature, "\nExpression")) +
            ggplot2::theme(
                legend.position = "right",
                legend.box = "vertical",
                legend.key = ggplot2::element_rect(fill = "white"),
                legend.background = ggplot2::element_blank()
            )

        # Lower facet - always scatter with expression coloring
        .expressionScatterFunc <- function(data, mapping, ...) {
            ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::geom_point(alpha = 0.7, size = 1.2,
                                    ggplot2::aes(color = .data[["Expression"]])) +
                ggplot2::scale_color_gradient(low = "gray", high = "blue") +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(
                        color = "black", fill = NA, linewidth = 0.5),
                    legend.position = "none")
        }

        # Blank facet function for diagonal and upper
        .blankFunc <- function(data, mapping, ...) {
            ggplot2::ggplot() +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    panel.border = ggplot2::element_rect(
                        color = "black", fill = NA, linewidth = 0.5),
                    legend.position = "none",
                    axis.text = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_blank(),
                    axis.title = ggplot2::element_blank(),
                    panel.grid = ggplot2::element_blank())
        }

        # Create pairs plot using GGally
        plot_obj <- suppressMessages(
            GGally::ggpairs(
                pc_df,
                columns = seq_len(length(pc_subset)),
                mapping = ggplot2::aes(color = .data[["Expression"]]),
                lower = list(continuous = .expressionScatterFunc),
                upper = list(continuous = .blankFunc),
                diag = list(continuous = .blankFunc),
                progress = FALSE,
                legend = GGally::grab_legend(legend_plot)
            )
        )

        # Add consistent theming
        plot_obj <- plot_obj +
            ggplot2::theme(
                strip.background = ggplot2::element_rect(
                    fill = "white", color = "black", linewidth = 0.5),
                strip.text = ggplot2::element_text(color = "black")
            )

        # Add title with cell type info if filtered
        if (!is.null(cell_type_col) && !is.null(cell_types)) {
            plot_obj <- plot_obj +
                ggplot2::labs(title = paste0(feature, " Expression in ",
                                             paste(cell_types, collapse = ", "),
                                             " Cells"))
        }
    }

    return(plot_obj)
}
