#' @title Visualization of gene sets or pathway scores on dimensional reduction plot
#'
#' @description
#' Plot gene sets or pathway scores on PCA, TSNE, or UMAP. Single cells are color-coded by scores of gene sets or pathways.
#'
#' @details
#' This function plots gene set scores on reduced dimensions such as PCA, t-SNE, or UMAP.
#' It extracts the reduced dimensions from the provided \code{\linkS4class{SingleCellExperiment}} object.
#' Gene set scores are visualized as a scatter plot with colors indicating the scores.
#' For PCA, the function automatically includes the percentage of variance explained
#' in the plot's legend.
#'
#' @param sce_object An object of class \code{\linkS4class{SingleCellExperiment}} containing numeric expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param cell_type_col The column name in the \code{colData} of \code{sce_object} that identifies the cell types.
#' @param method A character string indicating the method for visualization ("PCA", "TSNE", or "UMAP").
#' @param score_col A character string representing the name of the score_col (score) in the colData(sce_object) to plot.
#' @param pc_subset An optional vector specifying the principal components (PCs) to include in the plot if method = "PCA".
#'        Default is 1:5.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2000.
#'
#' @return A ggplot2 object representing the gene set scores plotted on the specified reduced dimensions.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("query_data")
#'
#' # Plot gene set scores on PCA
#' plotGeneSetScores(sce_object = query_data,
#'                   method = "PCA",
#'                   score_col = "gene_set_scores",
#'                   pc_subset = 1:5,
#'                   cell_types = "CD8",
#'                   cell_type_col = "SingleR_annotation")
#'
# Function to plot the gene set scores
plotGeneSetScores <- function(sce_object,
                              cell_type_col,
                              method = c("PCA", "TSNE", "UMAP"),
                              score_col,
                              pc_subset = 1:5,
                              cell_types = NULL,
                              max_cells = 2000) {

    # Check standard input arguments
    argumentCheck(query_data = sce_object,
                  pc_subset_query = pc_subset,
                  query_cell_type_col = cell_type_col,
                  max_cells_query = max_cells)

    # Convert cell type columns to character if needed
    sce_object <- convertColumnsToCharacter(sce_object = sce_object,
                                            convert_cols = cell_type_col)

    # Select cell types
    cell_types <- selectCellTypes(query_data = sce_object,
                                  query_cell_type_col = cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Match method argument
    method <- match.arg(method)

    # Store PCA attributes BEFORE downsampling (if method is PCA)
    pca_percent_var <- NULL
    if (method == "PCA" && "PCA" %in% reducedDimNames(sce_object)) {
        pca_attrs <- attributes(reducedDim(sce_object, "PCA"))
        if (!is.null(pca_attrs[["percentVar"]])) {
            pca_percent_var <- pca_attrs[["percentVar"]]
        }
    }

    # Downsample SCE object
    sce_object <- downsampleSCE(sce_object = sce_object,
                               max_cells = max_cells,
                               cell_types = cell_types,
                               cell_type_col = cell_type_col)

    # Check if dimension reduction method is present in reference's reducedDims
    if (!(method %in% names(reducedDims(sce_object)))) {
        stop("\'sce_object\' must have pre-computed \'",
             method,
             "\' in \'reducedDims\'.")
    }

    # Check if score_col is a valid column name in sce_object
    if (!score_col %in% colnames(colData(sce_object))) {
        stop("score_col: '", score_col,
             "' is not a valid column name in sce_object.")
    }

    # Handle cell type filtering
    filtered_to_specific_types <- FALSE
    if (!is.null(cell_type_col)) {
        # Check if cell type column exists
        if (!cell_type_col %in% colnames(colData(sce_object))) {
            stop("Specified cell_type_col does not exist in colData.")
        }

        # Get available cell types
        available_cell_types <- unique(colData(sce_object)[[cell_type_col]])

        # If specific cell types are requested, filter to those
        if (!is.null(cell_types)) {
            # Check if requested cell types exist
            missing_types <- setdiff(cell_types, available_cell_types)
            if (length(missing_types) > 0) {
                warning("The following cell types were not found: ",
                        paste(missing_types, collapse = ", "))
            }

            # Filter to cells of specified types
            cell_mask <- colData(sce_object)[[cell_type_col]] %in%
                cell_types
            filtered_to_specific_types <- TRUE
        } else {
            # Use all cells if no specific types requested
            cell_mask <- rep(TRUE, ncol(sce_object))
            cell_types <- available_cell_types
        }

        # Filter the object
        sce_object <- sce_object[, cell_mask]
    }

    # Extract scores (after potential filtering)
    scores <- colData(sce_object)[[score_col]]

    if(method %in% c("TSNE", "UMAP")){

        # Extract dimension reduction coordinates from SingleCellExperiment object
        reduction <- reducedDim(sce_object, method)

        # Prepare data for plotting
        df <- data.frame(Dim1 = reduction[, 1], Dim2 = reduction[, 2],
                         Scores = scores)

        # Add cell type information if available
        if (!is.null(cell_type_col)) {
            df[["CellType"]] <- colData(sce_object)[[cell_type_col]]
        }

        # Create the plot object with the original color scheme
        plot_obj <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Dim1"]],
                                                     y = .data[["Dim2"]])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["Scores"]]),
                                alpha = 0.7, size = 1.2) +
            ggplot2::scale_color_gradientn(
                colors = c("#2171B5", "#8AABC1", "#FFEDA0", "#E6550D"),
                values = seq(0, 1, by = 1/3),
                limits = c(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE)),
                name = "Scores") +
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
                ggplot2::ggtitle(paste0("Gene Set Scores in ",
                                        paste(cell_types, collapse = ", "),
                                        " Cells"))
        }

    } else if (method == "PCA"){

        # Check input for pc_subset
        if(!all(pc_subset %in% seq_len(ncol(reducedDim(sce_object, "PCA")))))
            stop("\'pc_subset\' is out of range.")

        # PCA data
        plot_mat <- reducedDim(sce_object, "PCA")[, pc_subset]

        # Create PC column names with variance explained (always show percentages if available)
        if (!is.null(pca_percent_var) && length(pca_percent_var) >= max(pc_subset)) {
            plot_names <- paste0(
                "PC", pc_subset, " (",
                sprintf("%.1f%%", pca_percent_var[pc_subset]), ")")
        } else {
            # Fallback if percentVar is not available
            plot_names <- paste0("PC", pc_subset)
        }

        # Create a new data frame with selected PCs
        pc_df <- data.frame(matrix(0, nrow = nrow(plot_mat),
                                   ncol = length(pc_subset)))
        colnames(pc_df) <- plot_names

        for (i in 1:length(pc_subset)) {
            pc_df[, i] <- plot_mat[, i]
        }

        # Add scores data
        pc_df[["Scores"]] <- scores

        # Add cell type information if available
        if (!is.null(cell_type_col)) {
            pc_df[["CellType"]] <- colData(sce_object)[[cell_type_col]]
        }

        # Create a simple plot to extract the legend with original color scheme
        legend_plot <- ggplot2::ggplot(pc_df,
                                       ggplot2::aes(x = pc_df[,1],
                                                    y = pc_df[,2])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["Scores"]])) +
            ggplot2::scale_color_gradientn(
                colors = c("#2171B5", "#8AABC1", "#FFEDA0", "#E6550D"),
                values = seq(0, 1, by = 1/3),
                limits = c(min(scores, na.rm = TRUE), max(scores,
                                                          na.rm = TRUE)),
                name = "Scores") +
            ggplot2::theme(
                legend.position = "right",
                legend.box = "vertical",
                legend.key = ggplot2::element_rect(fill = "white"),
                legend.background = ggplot2::element_blank()
            )

        # Lower facet - always scatter with score coloring
        .scoresScatterFunc <- function(data, mapping, ...) {
            ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::geom_point(alpha = 0.7, size = 1.2,
                                    ggplot2::aes(color = .data[["Scores"]])) +
                ggplot2::scale_color_gradientn(
                    colors = c("#2171B5", "#8AABC1", "#FFEDA0", "#E6550D"),
                    values = seq(0, 1, by = 1/3),
                    limits = c(min(scores, na.rm = TRUE), max(scores,
                                                              na.rm = TRUE))) +
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
                mapping = ggplot2::aes(color = .data[["Scores"]]),
                lower = list(continuous = .scoresScatterFunc),
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
                ggplot2::labs(title = paste0("Gene Set Scores in ",
                                             paste(cell_types,
                                                   collapse = ", "),
                                             " Cells"))
        }
    }

    return(plot_obj)
}
