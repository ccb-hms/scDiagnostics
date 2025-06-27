#' @title Plot Reference and Query Cell Types using MDS
#'
#' @description
#' This function facilitates the assessment of similarity between reference and query datasets
#' through Multidimensional Scaling (MDS) scatter plots. It allows the visualization of cell types,
#' color-coded with user-defined custom colors, based on a dissimilarity matrix computed from a
#' user-selected gene set.
#'
#' @details
#' To evaluate dataset similarity, the function selects specific subsets of cells from
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
                            assay_name = "logcounts") {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  assay_name = assay_name)

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Subset data
    query_data <- query_data[, which(
        query_data[[query_cell_type_col]] %in% cell_types)]
    reference_data <- reference_data[, which(
        reference_data[[ref_cell_type_col]] %in% cell_types)]

    # Extract assay matrices
    query_assay <- as.matrix(assay(query_data, assay_name))
    ref_assay <- as.matrix(assay(reference_data, assay_name))

    # Compute correlation and dissimilarity matrix
    df <- cbind(query_assay, ref_assay)
    corMat <- cor(df, method = "spearman")
    disMat <- (1 - corMat)
    cmd <- data.frame(cmdscale(disMat), c(rep("Query", ncol(query_assay)),
                                          rep("Reference", ncol(ref_assay))),
                      c(query_data[[query_cell_type_col]],
                        reference_data[[ref_cell_type_col]]))
    colnames(cmd) <- c("Dim1", "Dim2", "dataset", "cellType")
    cmd <- na.omit(cmd)
    cmd$cell_type_dataset <- paste(cmd[["dataset"]], cmd[["cellType"]], sep = " ")

    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)),
                                rep(sort(cell_types), each = 2))
    cmd$cell_type_dataset <- factor(cmd$cell_type_dataset,
                                    levels = order_combinations)

    # Define the colors for cell types
    cell_type_colors <- generateColors(order_combinations, paired = TRUE)

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
