#' @title Plot Cosine Similarities Between Cells and PCs
#'
#' @description The S3 plot method creates a heatmap plot to visualize the
#' cosine similarities between cells and principal components (PCs).
#'
#' @details The S3 plot method reshapes the input data frame to create a long
#' format suitable for plotting as a heatmap. It then creates a heatmap plot
#' using ggplot2, where the x-axis represents the PCs, the y-axis represents the
#' cells, and the color intensity represents the cosine similarity values.
#'
#' @param x An object of class `calculateCellSimilarityPCA` containing a
#' dataframe of cosine similarity values between cells and PCs.
#' @param pc_subset A numeric vector specifying the subset of principal
#' components to include in the plot. Default is 1:5.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @keywords internal
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the
#' cosine similarity heatmap.
#'
#' @export
#'
#' @author Anthony Christidis,
#' \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateCellSimilarityPCA}}
#'
#' @rdname calculateCellSimilarityPCA
#'
# Function to plot cosine similarities between cells and PCs
plot.calculateCellSimilarityPCAObject <- function(x,
                                                  pc_subset = 1:5,
                                                  ...) {
    # Subset data
    x <- x[, paste0("PC", pc_subset)]

    # Reshape the matrix into long format. Row-major flattening via
    # t()/as.vector() reproduces exactly the values and order the original
    # nested loop (outer over cells, inner over PCs) produced, but without
    # the O(n^2) cost of repeatedly growing vectors with c().
    if (is.matrix(x)) {
        n_cells <- nrow(x)
        n_pcs <- ncol(x)
        cell_names <- rep(rownames(x), each = n_pcs)
        pc_names <- rep(colnames(x), times = n_cells)
        cosine_values <- as.vector(t(x))
    } else {
        cell_names <- c()
        pc_names <- c()
        cosine_values <- c()
    }

    # Create a data frame with the reshaped data
    cosine_long <- data.frame(
        Cell = factor(cell_names,
            levels = rev(rownames(x))
        ),
        PC = pc_names, CosineSimilarity = cosine_values
    )

    # Create the heatmap plot
    plot <- ggplot2::ggplot(cosine_long, ggplot2::aes(
        x = .data[["PC"]],
        y = .data[["Cell"]],
        fill = .data[["CosineSimilarity"]]
    )) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%.2f", .data[["CosineSimilarity"]])),
            size = 3
        ) +
        ggplot2::scale_fill_gradient2(
            low = "blue", mid = "white",
            high = "red", midpoint = 0,
            limits = c(-1, 1), space = "Lab",
            name = "Cosine Similarity"
        ) +
        ggplot2::labs(title = "Cosine Similarity Heatmap", x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1
            ),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )
    return(plot)
}
