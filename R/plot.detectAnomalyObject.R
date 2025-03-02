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

    # Helper function to get combined ranges from reference and query data
    get_combined_ranges <- function(x, cell_type, pc_subset) {
        x_ranges <- list()
        ref_data <- x[[cell_type]]$reference_mat_subset[, pc_subset, drop = FALSE]
        query_data <- NULL
        if(!is.null(x[[cell_type]]$query_mat_subset)) {
            query_data <- x[[cell_type]]$query_mat_subset[, pc_subset, drop = FALSE]
        }
        for(i in seq_along(pc_subset)) {
            pc_values <- ref_data[, i]
            if(!is.null(query_data)) {
                pc_values <- c(pc_values, query_data[, i])
            }
            x_ranges[[i]] <- range(pc_values)
        }
        return(x_ranges)
    }

    # Get combined ranges for all PCs
    combined_ranges <- get_combined_ranges(x, cell_type, pc_subset)

    # Filter and prepare data based on data type
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

    # Add variance explained to column names
    colnames(data_subset) <- paste0(
        "PC", pc_subset,
        " (", sprintf("%.1f%%", x[[cell_type]]$var_explained[pc_subset]), ")")

    # Create all possible PC pairs
    pc_names <- colnames(data_subset)
    pairs <- expand.grid(x = pc_names, y = pc_names)
    pairs <- pairs[pairs[["x"]] != pairs[["y"]], ]

    # Create data frame with PC pairs
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

    # Remove redundant pairs to avoid duplicate plots
    data_pairs <- data_pairs[as.numeric(data_pairs[["x"]]) <
                                 as.numeric(data_pairs[["y"]]),]

    # Add anomaly information
    data_pairs[["anomaly"]] <- factor(
        rep(anomaly, choose(length(pc_subset), 2)),
        levels = c("TRUE", "FALSE"))

    # Create background data using combined ranges
    pc_labels <- unique(c(data_pairs[["x"]], data_pairs[["y"]]))
    pc_combns <- combn(pc_subset, 2, simplify = FALSE)
    background_data <- data.frame("x_value" = numeric(0), "y_value" = numeric(0),
                                  "x" = character(0), "y" = character(0),
                                  anomaly_score = numeric(0))

    # Generate background grid for each PC combination
    for(i in seq_along(pc_combns)){
        pc_combn <- pc_combns[[i]]

        # Get pre-calculated ranges and add buffer
        x_range <- combined_ranges[[which(pc_subset == pc_combn[1])]]
        y_range <- combined_ranges[[which(pc_subset == pc_combn[2])]]
        x_buffer <- 0.05 * diff(x_range)
        y_buffer <- 0.05 * diff(y_range)

        # Create sequences for grid
        x_seq <- seq(x_range[1] - x_buffer, x_range[2] + x_buffer, length.out = 200)
        y_seq <- seq(y_range[1] - y_buffer, y_range[2] + y_buffer, length.out = 200)
        grid <- expand.grid(x = x_seq, y = y_seq)

        # Calculate anomaly scores
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

    # Calculate width and height for tiles based on data ranges with expanded limits
    tile_sizes <- lapply(pc_combns, function(pc_combn) {
        x_data <- background_data[background_data$x == pc_labels[which(pc_subset == pc_combn[1])], "x_value"]
        y_data <- background_data[background_data$y == pc_labels[which(pc_subset == pc_combn[2])], "y_value"]

        # Expand ranges by 5%
        x_range <- diff(range(x_data))
        y_range <- diff(range(y_data))
        x_buffer <- 0.05 * x_range
        y_buffer <- 0.05 * y_range

        data.frame(
            x = pc_labels[which(pc_subset == pc_combn[1])],
            y = pc_labels[which(pc_subset == pc_combn[2])],
            width = x_range/50,
            height = y_range/50
        )
    })
    tile_sizes <- do.call(rbind, tile_sizes)

    # Merge tile sizes with background data
    background_data <- merge(background_data, tile_sizes, by = c("x", "y"))

    # Create final plot
    anomaly_plot <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = background_data,
                           ggplot2::aes(x = .data[["x_value"]],
                                        y = .data[["y_value"]],
                                        fill = .data[["anomaly_score"]],
                                        width = .data[["width"]],
                                        height = .data[["height"]])) +
        ggplot2::scale_fill_gradient2(
            low = "lightgreen",
            mid = "lightyellow",
            high = "lightpink",
            midpoint = 0.5,
            guide = "none") +
        ggplot2::facet_grid(rows = ggplot2::vars(.data[["y"]]),
                            cols = ggplot2::vars(.data[["x"]]),
                            scales = "free") +
        ggplot2::geom_point(data = data_pairs,
                            ggplot2::aes(x = .data[["x_value"]],
                                         y = .data[["y_value"]],
                                         color = factor(.data[["anomaly"]])),
                            shape = 16,
                            size = 1.5,
                            alpha = 0.45) +
        ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                                    breaks = c("TRUE", "FALSE"),
                                    labels = c("TRUE", "FALSE"),
                                    name = "Anomalous") +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
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
                                     cell_type))

    return(anomaly_plot)
}
