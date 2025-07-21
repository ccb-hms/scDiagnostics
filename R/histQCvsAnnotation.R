#' @title Histograms: QC Stats and Annotation Scores Visualization
#'
#' @description
#' This function generates histograms for visualizing the distribution of quality control (QC) statistics and
#' annotation scores associated with cell types in single-cell genomic data.
#'
#' @details The particularly useful in the analysis of data from single-cell experiments,
#' where understanding the distribution of these metrics is crucial for quality assessment and
#' interpretation of cell type annotations.
#'
#' @param se_object  A \code{\linkS4class{SingleCellExperiment}} containing the single-cell
#' expression data and metadata.
#' @param cell_type_col The column name in the \code{colData} of \code{se_object}
#' that contains the cell type labels.
#' @param cell_types A vector of cell types to plot (e.g., c("T-cell", "B-cell")).
#' Defaults to \code{NULL}, which will include all the cells.
#' @param qc_col A column name in the \code{colData} of \code{se_object} that
#' contains the QC stats of interest.
#' @param score_col The column name in the \code{colData} of \code{se_object} that
#' contains the cell type scores.
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A object containing two histograms displayed side by side.
#' The first histogram represents the distribution of QC stats,
#' and the second histogram represents the distribution of annotation scores.
#'
#' @examples
# Load data
#' data("query_data")
#'
#' # Generate histograms
#' histQCvsAnnotation(se_object = query_data,
#'                    cell_type_col = "SingleR_annotation",
#'                    cell_types = c("CD4", "CD8"),
#'                    qc_col = "percent_mito",
#'                    score_col = "annotation_scores")
#'
#' histQCvsAnnotation(se_object = query_data,
#'                    cell_type_col = "SingleR_annotation",
#'                    cell_types = NULL,
#'                    qc_col = "percent_mito",
#'                    score_col = "annotation_scores")
#'
#' @export
#'
# Function to plot histogram of QC scores and annotation scores
histQCvsAnnotation <- function(se_object,
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
        stop("qc_col: '",
             qc_col,
             "' is not a valid column name in se_object.")
    }

    # Check if score_col is a valid column name in se_object
    if (!score_col %in% names(colData(se_object))) {
        stop("score_col: '",
             score_col,
             "' is not a valid column name in se_object.")
    }

    # Filter cells based on cell_types if specified
    if (!is.null(cell_types)) {
        index <- which(se_object[[cell_type_col]] %in% cell_types)
        se_object <- se_object[, index]
    }

    # Extract QC stats, scores, and cell_types
    qc_stats <- se_object[[qc_col]]
    cell_type_scores <- se_object[[score_col]]

    # Combine QC stats, scores, and cell_types into a data frame
    data <- data.frame(stats = c(qc_stats, cell_type_scores),
                       Metric = c(rep("QC Statistics", length(qc_stats)),
                                  rep("Annotation Scores",
                                      length(cell_type_scores))))
    data[["Metric"]] <- factor(data[["Metric"]],
                               levels = c("QC Statistics", "Annotation Scores"))

    # Create histogram plots
    hist_plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["stats"]],
                                                    fill = .data[["Metric"]])) +
        ggplot2::geom_histogram(bins = 30, alpha = 0.5, position = "identity",
                                color = "gray50") +
        ggplot2::facet_wrap(~ .data[["Metric"]], scales = "free") +
        ggplot2::labs(x = "", y = "Frequency", fill = "Metric") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "white", color = "black", linewidth = 0.5),
            legend.position = "none",
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold",
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10))

    # Return the list of plots
    return(hist_plot)
}
