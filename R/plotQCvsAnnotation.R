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
#' @param se_object A \code{\linkS4class{SingleCellExperiment}} containing the single-cell
#' expression data and metadata.
#' @param cell_type_col The column name in the \code{colData} of \code{se_object}
#' that contains the cell type labels.
#' @param cell_types A vector of cell type labels to plot (e.g., c("T-cell", "B-cell")).
#' Defaults to \code{NULL}, which will include all the cells.
#' @param qc_col A column name in the \code{colData} of \code{se_object} that
#' contains the QC stats of interest.
#' @param score_col The column name in the \code{colData} of \code{se_object} that
#' contains the cell type annotation scores.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
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
#' p1 <- plotQCvsAnnotation(se_object = qc_data_subset,
#'                          cell_type_col = "SingleR_annotation",
#'                          cell_types = NULL,
#'                          qc_col = "total",
#'                          score_col = "annotation_scores")
#'p1 + ggplot2::xlab("Library Size")
#'
#' @export
#'
# Function to plot QC score against annotation
plotQCvsAnnotation <- function(se_object,
                               cell_type_col,
                               cell_types = NULL,
                               qc_col,
                               score_col,
                               max_cells = 2500) {

    # Check standard input arguments
    argumentCheck(query_data = se_object,
                  query_cell_type_col = cell_type_col,
                  cell_types = cell_types)

    # Downsample SCE object
    se_object <- downsampleSCE(sce = se_object,
                               max_cells = max_cells)

    # Check if qc_col is a valid column name in se_object
    if (!qc_col %in% names(colData(se_object))) {
        stop("qc_col: '", qc_col, "' is not a valid column name in se_object.")
    }

    # Check if score_col is a valid column name in se_object
    if (!score_col %in% names(colData(se_object))) {
        stop("score_col: '", score_col, "' is not a valid column name in se_object.")
    }

    # Filter cells based on cell_types if specified
    if (!is.null(cell_types)) {
        se_object <- se_object[, which(se_object[[cell_type_col]] %in%
                                           cell_types)]
    }

    # Extract QC stats, scores, and cell_typess
    qc_stats <- se_object[[qc_col]]
    cell_types_scores <- se_object[[score_col]]
    cell_labels <- se_object[[cell_type_col]]

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
