#' @title Create Faceted Scatter Plots for Specified PC Combinations From \code{detectAnomaly} Object
#'
#' @description
#' This S3 plot method generates faceted scatter plots for specified principal component (PC) combinations
#' within an anomaly detection object. It visualizes the relationship between specified PCs, highlights
#' anomalies detected by the Isolation Forest algorithm, and provides a background gradient representing
#' anomaly scores.
#'
#' @details
#' The function extracts the specified PCs from the given anomaly detection object and generates
#' scatter plots for each pair of PCs. It uses \code{ggplot2} to create a faceted plot where each facet
#' represents a pair of PCs. The plot includes the following elements:
#'
#' 1. A background gradient: This represents the anomaly scores across the PC space. The gradient
#'    ranges from green (low anomaly scores) through yellow to red (high anomaly scores).
#' 2. Data points: Plotted over the gradient, with different shapes for normal and anomalous points.
#'    Anomalous points are represented as 'X', while normal points are solid circles.
#'
#' The background gradient provides a visual representation of the anomaly landscape, allowing for
#' intuitive interpretation of regions with high or low anomaly scores.
#'
#' @param x A list object containing the anomaly detection results from the \code{detectAnomaly} function.
#' Each element of the list should correspond to a cell type and contain \code{reference_mat_subset}, \code{query_mat_subset},
#' \code{var_explained}, and \code{anomaly}.
#' @param cell_type A character string specifying the cell type for which the plots should be generated. This should
#' be a name present in \code{x}. If NULL, the "Combined" cell type will be plotted. Default is NULL.
#' @param pc_subset A numeric vector specifying the indices of the PCs to be included in the plots. If NULL, all PCs
#' in \code{reference_mat_subset} will be included.
#' @param data_type A character string specifying whether to plot the "query" data or the "reference" data. Default is "query".
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 500
#' @param ... Additional arguments passed to the `isolation.forest` function.
#'
#' @keywords internal
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the PCA plots with anomalies highlighted.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{detectAnomaly}}
#'
#' @rdname detectAnomaly
#'
# Function to create faceted scatter plots for specified PC combinations
plot.detectAnomalyObject <- function(x,
                                     cell_type = NULL,
                                     pc_subset = NULL,
                                     data_type = c("query", "reference"),
                                     n_tree = 500,
                                     ...) {

    # Check if PCA was used for computations
    if(!("var_explained" %in% names(x[[names(x)[1]]])))
        stop("The plot function can only be used if \'n_components\' is not NULL.")

    # Check input for cell type
    if(is.null(cell_type)){
        cell_type <- "Combined"
    } else{
        if(!(cell_type %in% names(x)))
            stop("\'cell_type\' is not available in \'x\'.")
    }

    # Check input for pc_subset
    if(!is.null(pc_subset)){
        if(!all(pc_subset %in% seq_len(
            ncol(x[[cell_type]][["reference_mat_subset"]]))))
            stop("\'pc_subset\' is out of range.")
    } else{
        pc_subset <- seq_len(ncol(x[[cell_type]][["reference_mat_subset"]]))
    }

    # Check input for data_type
    data_type <- match.arg(data_type)

    # Filter data to include only specified PCs
    if(is.null(x[[cell_type]]$query_mat_subset) && data_type == "query"){
        stop("There is no query data available in the \'detectAnomaly\' object.")
    } else{
        if(data_type == "query"){
            data_subset <- x[[cell_type]]$query_mat_subset[, pc_subset,
                                                           drop = FALSE]
            anomaly <- x[[cell_type]]$query_anomaly

        } else if(data_type == "reference"){
            data_subset <- x[[cell_type]]$reference_mat_subset[, pc_subset,
                                                               drop = FALSE]
            anomaly <- x[[cell_type]]$reference_anomaly
        }
    }

    # Modify column names to include percentage of variance explained
    colnames(data_subset) <- paste0(
        "PC", pc_subset,
        " (", sprintf("%.1f%%", x[[cell_type]]$var_explained[pc_subset]), ")")

    # Create all possible pairs of specified PCs
    pc_names <- colnames(data_subset)
    pairs <- expand.grid(x = pc_names, y = pc_names)
    pairs <- pairs[pairs[["x"]] != pairs[["y"]], ]

    # Create a new data frame with all possible pairs of specified PCs
    data_pairs_list <- lapply(seq_len(nrow(pairs)), function(i) {
        x_col <- pairs[["x"]][i]
        y_col <- pairs[["y"]][i]
        data_frame <- data.frame(data_subset[, c(x_col, y_col)])
        colnames(data_frame) <- c("x_value", "y_value")
        data_frame[["x"]] <- x_col
        data_frame[["y"]] <- y_col
        return(data_frame)
    })
    data_pairs <- do.call(rbind, data_pairs_list)

    # Remove redundant data (to avoid duplicated plots)
    data_pairs <- data_pairs[as.numeric(data_pairs[["x"]]) <
                                 as.numeric(data_pairs[["y"]]),]

    # Add anomalies vector to data_pairs dataframe
    data_pairs[["anomaly"]] <- factor(
        rep(anomaly, choose(length(pc_subset), 2)),
        levels = c("TRUE", "FALSE"))

    # Create background data
    pc_labels <- unique(c(data_pairs[["x"]], data_pairs[["y"]]))
    pc_combns <- combn(pc_subset, 2, simplify = FALSE)
    background_data <- data.frame("x_value" = numeric(0), "y_value" = numeric(0),
                                  "x" = character(0), "y" = character(0),
                                  anomaly_score = numeric(0))

    for(pc_combn in pc_combns){

        # Range for PC values
        x_range <- range(
            data_pairs[data_pairs[["x"]] ==
                           pc_labels[pc_combn[1]], ][["x_value"]])
        y_range <- range(
            data_pairs[data_pairs[["y"]] ==
                           pc_labels[pc_combn[2]], ][["y_value"]])

        # Create grid for background
        x_seq <- seq(x_range[1] - 1, x_range[2] + 1, length.out = 200)
        y_seq <- seq(y_range[1] - 1, y_range[2] + 1, length.out = 200)
        grid <- expand.grid(x = x_seq, y = y_seq)

        # Adding isolation forest prediction
        isolation_forest <- isotree::isolation.forest(
            x[[cell_type]][["reference_mat_subset"]][, paste0("PC", pc_combn)],
            ntree = n_tree)
        grid <- cbind(grid,
                      pc_labels[pc_combn[1]],
                      pc_labels[pc_combn[2]],
                      predict(isolation_forest,
                              newdata = grid,
                              type = "score"))
        colnames(grid) <- c("x_value", "y_value", "x", "y", "anomaly_score")
        background_data <- rbind(background_data, grid)
    }

    # Create the ggplot object with facets
    anomaly_plot <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = background_data,
                           ggplot2::aes(x = .data[["x_value"]],
                                        y = .data[["y_value"]],
                                        fill = .data[["anomaly_score"]]),
                           width = 0.15,
                           height = 0.15) +
        ggplot2::scale_fill_gradient2(
            low = "lightgreen", mid = "lightyellow", high = "lightpink",
            midpoint = mean(range(background_data[["anomaly_score"]])),
            guide = "none") +
        ggplot2::facet_grid(rows = ggplot2::vars(.data[["y"]]),
                            cols = ggplot2::vars(.data[["x"]]),
                            scales = "free") +
        ggplot2::geom_point(data = data_pairs, ggplot2::aes(
            x = .data[["x_value"]], y = .data[["y_value"]],
            color = factor(.data[["anomaly"]])),
            shape = 16,
            size = 1.5,
            alpha = 0.45) +
        ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                                    breaks = c("TRUE", "FALSE"),
                                    labels = c("TRUE", "FALSE"),
                                    name = "Anomalous") +
        ggplot2::xlab("") + ggplot2::ylab("") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = 14,
                                                          face = "bold",
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 12),
                       axis.text = ggplot2::element_text(size = 10),
                       legend.position = "right") +
        ggplot2::labs(title = paste0("Isolation Forest Anomaly Plot: ",
                                     cell_type))  +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(expand = FALSE)
    return(anomaly_plot)
}
