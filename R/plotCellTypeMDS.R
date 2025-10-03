#' @title Plot Reference and Query Cell Types using MDS
#'
#' @description
#' This function facilitates the assessment of similarity between reference and query datasets
#' through Multidimensional Scaling (MDS) scatter plots. It allows the visualization of cell types,
#' color-coded with user-defined custom colors, based on a dissimilarity matrix computed from a
#' user-selected gene set. If MDS coordinates are precomputed in reducedDims, they will be used;
#' otherwise, MDS will be computed from scratch.
#'
#' @details
#' The function first checks if MDS coordinates are available in the reducedDims of both datasets.
#' If precomputed MDS is found, it uses those coordinates directly for visualization.
#'
#' If MDS is not precomputed, the function selects specific subsets of cells from
#' both reference and query datasets. It then calculates Spearman correlations between gene expression profiles,
#' deriving a dissimilarity matrix. This matrix undergoes Classical Multidimensional Scaling (MDS) for
#' visualization, presenting cell types in a scatter plot, distinguished by colors defined by the user.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} containing the single-cell
#' expression data and metadata.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing the single-cell
#' expression data and metadata.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#'
#' @return A ggplot object representing the MDS scatter plot with cell type coloring.
#'
#' @references
#' \itemize{
#' \item Kruskal, J. B. (1964). "Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis". *Psychometrika*, 29(1), 1-27. doi:10.1007/BF02289565.
#' \item Borg, I., & Groenen, P. J. F. (2005). *Modern multidimensional scaling: Theory and applications* (2nd ed.). Springer Science & Business Media. doi:10.1007/978-0-387-25975-1.
#' }
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Generate the MDS scatter plot with cell type coloring
#' mds_plot <- plotCellTypeMDS(query_data = query_data,
#'                             reference_data = reference_data,
#'                             cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid")[1:4],
#'                             query_cell_type_col = "SingleR_annotation",
#'                             ref_cell_type_col = "expert_annotation")
#' mds_plot
#'
#' @importFrom stats cmdscale cor
#' @importFrom SummarizedExperiment assay
#'
# Function to visualize cell types in MDS plot
plotCellTypeMDS <- function(query_data,
                            reference_data,
                            query_cell_type_col,
                            ref_cell_type_col,
                            cell_types = NULL,
                            assay_name = "logcounts",
                            max_cells_query = 5000,
                            max_cells_ref = 5000) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cell_type_col)
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = 10)

    # Check if MDS is precomputed in both datasets
    query_has_mds <- "MDS" %in% reducedDimNames(query_data)
    ref_has_mds <- "MDS" %in% reducedDimNames(reference_data)
    use_precomputed_mds <- query_has_mds && ref_has_mds

    if (use_precomputed_mds) {
        message("Using precomputed MDS coordinates from reducedDims.")

        # Downsample using precomputed MDS (coordinates will be preserved through subsetting)
        query_data <- downsampleSCE(sce_object = query_data,
                                    max_cells = max_cells_query,
                                    cell_types = cell_types,
                                    cell_type_col = query_cell_type_col)
        reference_data <- downsampleSCE(sce_object = reference_data,
                                        max_cells = max_cells_ref,
                                        cell_types = cell_types,
                                        cell_type_col = ref_cell_type_col)

        # Extract precomputed MDS coordinates
        query_mds <- reducedDim(query_data, "MDS")
        ref_mds <- reducedDim(reference_data, "MDS")

        # Validate MDS dimensions (ensure at least 2D)
        if (ncol(query_mds) < 2 || ncol(ref_mds) < 2) {
            warning("Precomputed MDS has fewer than 2 dimensions. Computing MDS from scratch.")
            use_precomputed_mds <- FALSE
        }

        if (use_precomputed_mds) {
            # Create data frame from precomputed coordinates
            cmd <- data.frame(
                Dim1 = c(ref_mds[, 1], query_mds[, 1]),
                Dim2 = c(ref_mds[, 2], query_mds[, 2]),
                dataset = c(rep("Reference", nrow(ref_mds)), rep("Query", nrow(query_mds))),
                cellType = c(reference_data[[ref_cell_type_col]], query_data[[query_cell_type_col]]),
                stringsAsFactors = FALSE
            )
        }
    }

    # Fall back to computing MDS if not precomputed or validation failed
    if (!use_precomputed_mds) {
        if (query_has_mds || ref_has_mds) {
            message("MDS not available in both datasets or validation failed. Computing MDS from expression data.")
        } else {
            message("Computing MDS from expression data.")
        }

        # Downsample query and reference data (with cell type filtering)
        query_data <- downsampleSCE(sce_object = query_data,
                                    max_cells = max_cells_query,
                                    cell_types = cell_types,
                                    cell_type_col = query_cell_type_col)
        reference_data <- downsampleSCE(sce_object = reference_data,
                                        max_cells = max_cells_ref,
                                        cell_types = cell_types,
                                        cell_type_col = ref_cell_type_col)

        # Extract assay matrices
        query_assay <- as.matrix(assay(query_data, assay_name))
        ref_assay <- as.matrix(assay(reference_data, assay_name))

        # Compute correlation and dissimilarity matrix
        df <- cbind(query_assay, ref_assay)
        corMat <- cor(df, method = "spearman")
        disMat <- (1 - corMat)

        # Compute MDS coordinates
        mds_coords <- cmdscale(disMat)

        # Create data frame
        cmd <- data.frame(
            Dim1 = mds_coords[, 1],
            Dim2 = mds_coords[, 2],
            dataset = c(rep("Query", ncol(query_assay)), rep("Reference", ncol(ref_assay))),
            cellType = c(query_data[[query_cell_type_col]], reference_data[[ref_cell_type_col]]),
            stringsAsFactors = FALSE
        )
    }

    # Clean up data and create combined cell type-dataset labels
    cmd <- na.omit(cmd)
    cmd[["cell_type_dataset"]] <- paste(cmd[["dataset"]], cmd[["cellType"]], sep = " ")

    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)),
                                rep(sort(cell_types), each = 2))
    cmd[["cell_type_dataset"]] <- factor(cmd[["cell_type_dataset"]],
                                         levels = order_combinations)

    # Define the colors for cell types
    cell_type_colors <- generateColors(order_combinations, paired = TRUE)

    # Create the plot
    mds_plot <- ggplot2::ggplot(cmd,
                                ggplot2::aes(
                                    x = .data[["Dim1"]],
                                    y = .data[["Dim2"]],
                                    color = .data[["cell_type_dataset"]])) +
        ggplot2::geom_point(alpha = 0.5, size = 1) +
        ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Types") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14, face = "bold",
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10)) +
        ggplot2::guides(color = ggplot2::guide_legend(title = "Cell Types"))

    return(mds_plot)
}
