#' Scatter plot: QC stats vs Cell Type Annotation Scores
#'
#' Creates a scatter plot to visualize the relationship between QC stats (e.g., library size)
#' and cell type annotation scores for one or more cell types.
#'
#' @details This function generates a scatter plot to explore the relationship between various quality
#' control (QC) statistics, such as library size and mitochondrial percentage, and cell type
#' annotation scores. By examining these relationships, users can assess whether specific QC
#' metrics, systematically influence the confidence in cell type annotations,
#' which is essential for ensuring reliable cell type annotation.
#'
#' @param sce_object A \code{\linkS4class{SingleCellExperiment}} containing the single-cell
#' expression data and metadata.
#' @param cell_type_col The column name in the \code{colData} of \code{sce_object}
#' that contains the cell type labels.
#' @param cell_types A vector of cell type labels to plot (e.g., c("T-cell", "B-cell")).
#' Defaults to \code{NULL}, which will include all the cells.
#' @param qc_col A column name in the \code{colData} of \code{sce_object} that
#' contains the QC stats of interest.
#' @param score_col The column name in the \code{colData} of \code{sce_object} that
#' contains the cell type annotation scores.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 5000.
#'
#' @return A ggplot object displaying a scatter plot of QC stats vs annotation scores,
#'         where each point represents a cell, color-coded by its cell type.
#'
#' @examples
#' # Load data
#' data("qc_data")
#'
#' # Remove cell types with very few cells
#' qc_data_subset <- qc_data[, !(qc_data$SingleR_annotation
#'                               %in% c("Chondrocytes", "DC",
#'                                      "Neurons","Platelets"))]
#'
# Create a scatter plot between library size and annotation scores
#' p1 <- plotQCvsAnnotation(sce_object = qc_data_subset,
#'                          cell_type_col = "SingleR_annotation",
#'                          cell_types = NULL,
#'                          qc_col = "total",
#'                          score_col = "annotation_scores")
#'p1 + ggplot2::xlab("Library Size")
#'
#' @export
#'
# Function to plot QC score against annotation
plotQCvsAnnotation <- function(sce_object,
                               cell_type_col,
                               cell_types = NULL,
                               qc_col,
                               score_col,
                               max_cells = 5000) {

    # Check standard input arguments
    argumentCheck(query_data = sce_object,
                  query_cell_type_col = cell_type_col,
                  max_cells_query = max_cells)

    # Convert cell type columns to character if needed
    sce_object <- convertColumnsToCharacter(sce_object = sce_object,
                                            convert_cols = cell_type_col)

    # Downsample SCE object
    sce_object <- downsampleSCE(sce_object = sce_object,
                                max_cells = max_cells,
                                cell_types = cell_types,
                                cell_type_col = cell_type_col)

    # Check if qc_col is a valid column name in sce_object
    if (!qc_col %in% names(colData(sce_object))) {
        stop("qc_col: '", qc_col, "' is not a valid column name in sce_object.")
    }

    # Check if score_col is a valid column name in sce_object
    if (!score_col %in% names(colData(sce_object))) {
        stop("score_col: '", score_col, "' is not a valid column name in sce_object.")
    }

    # Select cell types
    cell_types <- selectCellTypes(query_data = sce_object,
                                  reference_data = NULL,
                                  query_cell_type_col = cell_type_col,
                                  ref_cell_type_col = NULL,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = 10)

    # Extract QC stats, scores, and cell_typess
    qc_stats <- sce_object[[qc_col]]
    cell_types_scores <- sce_object[[score_col]]
    cell_labels <- sce_object[[cell_type_col]]

    # Combine QC stats, scores, and labels into a data frame
    data <- data.frame(QCStats = qc_stats,
                       Scores = cell_types_scores,
                       cell_type = cell_labels)

    # Define the colors for cell types
    cell_type_colors <- generateColors(sort(unique(cell_labels)),
                                       paired = FALSE)

    # Create a scatter plot with color-coded points based on cell types or labels
    qc_plot <- ggplot2::ggplot(data,
                               ggplot2::aes(x = .data[["QCStats"]],
                                            y = .data[["Scores"]],
                                            color = .data[["cell_type"]])) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = cell_type_colors,
                                    name = "Cell Type") +
        ggplot2::xlab("QC stats") +
        ggplot2::ylab("Annotation Scores") +
        ggplot2::labs(color = "Cell Type") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold", hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10))

    return(qc_plot)
}
