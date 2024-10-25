#' @title Plot Density of Wasserstein Distances for Null Distribution
#'
#' @description
#' This function generates a density plot of Wasserstein distances for the null distribution
#' of a `calculateWassersteinDistanceObject`. Additionally, it overlays lines representing
#' the significance threshold and the reference-query distance.
#'
#' @details
#' The density plot visualizes the distribution of Wasserstein distances calculated among
#' reference samples, representing the null distribution. A vertical line marks the
#' significance threshold based on the specified \code{alpha}. Another line indicates the
#' mean Wasserstein distance between the reference and query datasets.
#'
#' @param x A list object containing the Wasserstein distance results from the \code{calculateWassersteinDistance} function.
#' @param alpha A numeric value specifying the significance level for thresholding. Default is 0.05.
#' @param ... Additional arguments for future extensions.
#'
#' @keywords internal
#'
#' @return A ggplot2 object representing the ridge plots of Wasserstein distances with annotated p-value.
#'
#' @references
#' Schuhmacher, D., Bernhard, S., & Book, M. (2019). "A Review of Approximate Transport in Machine Learning".
#' In *Journal of Machine Learning Research* (Vol. 20, No. 117, pp. 1-61).
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{calculateWassersteinDistance}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Extract CD4 cells
#' ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
#' query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]
#'
#' # Selecting highly variable genes (can be customized by the user)
#' ref_top_genes <- scran::getTopHVGs(ref_data_subset, n = 500)
#' query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_top_genes, query_top_genes)
#' ref_data_subset <- ref_data_subset[common_genes,]
#' query_data_subset <- query_data_subset[common_genes,]
#'
#' # Run PCA on reference data
#' ref_data_subset <- scater::runPCA(ref_data_subset)
#'
#' # Compute Wasserstein null distribution using reference data and observed distances with query data
#' wasserstein_data <- calculateWassersteinDistance(query_data = query_data_subset,
#'                                                  reference_data = ref_data_subset,
#'                                                  query_cell_type_col = "expert_annotation",
#'                                                  ref_cell_type_col = "expert_annotation",
#'                                                  pc_subset = 1:5,
#'                                                  n_resamples = 100)
#' plot(wasserstein_data)
#'
# Function to generate density of Wasserstein distances under null distribution
plot.calculateWassersteinDistanceObject <- function(
        x,
        alpha = 0.05,
        ...){

    # Input check for alpha
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
        stop("\'alpha\' must be a positive number greater than 0 and less than 1.")
    }

    # Visualize results
    threshold_text <- bquote(paste(
        "Signifiance Threshold (", alpha, " = ", .(alpha), ")"))
    vline_data <- data.frame(xintercept = c(quantile(x$null_dist, 1 - alpha),
                                            x$query_dist),
                             line_type = c("Signifiance Threshold",
                                           "Reference-Query Distance"))
    density_plot <- ggplot2::ggplot(data.frame(x$null_dist),
                                    ggplot2::aes(x = x$null_dist)) +
        ggplot2::geom_density(alpha = 0.7, fill = "#00BBC4") +
        ggplot2::labs(title = paste0(
            "Density of Wasserstein Distances For Reference Distribution of ",
            x$cell_type),
            x = "Wasserstein Distances", y = "Density") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "right",
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14, face = "bold",
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10)) +
        ggplot2::geom_vline(data = vline_data,
                            ggplot2::aes(xintercept = .data[["xintercept"]],
                                         linetype = .data[["line_type"]]),
                            color = "black", linewidth = c(1, 1)) +
        ggplot2::scale_linetype_manual(
            name = NULL,
            values = c("Signifiance Threshold" = "solid",
                       "Reference-Query Distance" = "dashed"),
            labels = c("Reference-Query Distance", threshold_text)) +
        ggplot2::guides(linetype = ggplot2::guide_legend(
            nrow = 2, override.aes = list(color = "black", size = 0.5),
            direction = "horizontal",
            keywidth = ggplot2::unit(1, "line"),
            keyheight = ggplot2::unit(1.5, "line")))
    return(density_plot)
}

