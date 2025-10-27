#' @title Plot Heatmap of Similarities Between Principal Components
#'
#' @description
#' The S3 plot method generates a heatmap to visualize the similarities between
#' principal components from the output of the \code{comparePCA} function.
#'
#' @details
#' The S3 plot method creates an enhanced heatmap visualization with options to display
#' statistical significance and similarity values. The heatmap uses a blue-white-red
#' color gradient for similarity values, and optionally overlays significance indicators.
#'
#' @param x A \code{comparePCAObject} output from the \code{comparePCA} function.
#' @param show_values Logical, whether to display similarity values on the heatmap. Default is TRUE.
#' @param show_significance Logical, whether to display significance indicators (requires permutation test). Default is TRUE.
#' @param significance_threshold Numeric, p-value threshold for significance. Default is 0.05.
#' @param color_limits Numeric vector of length 2 specifying color scale limits. If NULL, uses data range.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A \code{ggplot} object representing the heatmap of similarities.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{comparePCA}}
#'
#' @rdname comparePCA
#'
# Function to plot the comparison of the PCs of the SCE objects
plot.comparePCAObject <- function(x,
                                  show_values = TRUE,
                                  show_significance = TRUE,
                                  significance_threshold = 0.05,
                                  color_limits = NULL,
                                  ...) {

    # Extract similarity matrix from the object
    if (is.list(x) && "similarity_matrix" %in% names(x)) {
        similarity_matrix <- x[["similarity_matrix"]]
        p_values <- x[["p_values"]]
        metric_name <- x[["metric"]]
    } else {
        # Backward compatibility with old format
        similarity_matrix <- x
        p_values <- NULL
        metric_name <- "cosine"
    }

    # Convert the matrix to a data frame
    similarity_df <- data.frame(
        Ref = factor(rep(rownames(similarity_matrix),
                         each = ncol(similarity_matrix)),
                     levels = rev(rownames(similarity_matrix))),
        Query = factor(rep(colnames(similarity_matrix),
                           times = nrow(similarity_matrix)),
                       levels = colnames(similarity_matrix)),
        Similarity = as.vector(similarity_matrix)
    )

    # Add significance information if available
    if (!is.null(p_values) && show_significance) {
        similarity_df[["P_Value"]] <- as.vector(p_values)
        similarity_df[["Significant"]] <- similarity_df[["P_Value"]] <
            significance_threshold
        similarity_df[["Significance_Symbol"]] <-
            ifelse(similarity_df[["Significant"]], "*", "")

        # Create combined labels with values and significance BEFORE ggplot
        if (show_values) {
            similarity_df[["Label"]] <- paste0(
                sprintf("%.2f", similarity_df[["Similarity"]]),
                similarity_df[["Significance_Symbol"]]
            )
        }
    } else {
        similarity_df[["Significant"]] <- FALSE
        similarity_df[["Significance_Symbol"]] <- ""

        # Create simple labels if just showing values
        if (show_values) {
            similarity_df[["Label"]] <-
                sprintf("%.2f", similarity_df[["Similarity"]])
        }
    }

    # Determine color limits
    if (is.null(color_limits)) {
        max_abs <- max(abs(range(similarity_matrix, na.rm = TRUE)))
        color_limits <- c(-max_abs, max_abs)
    }

    # Create appropriate title based on metric
    title_text <- switch(metric_name,
                         "cosine" = "Heatmap of Cosine Similarities Between PCs",
                         "correlation" = "Heatmap of Correlations Between PCs",
                         "Heatmap of Similarities Between PCs")

    # Create the base heatmap
    pc_plot <- ggplot2::ggplot(
        similarity_df,
        ggplot2::aes(x = .data[["Query"]],
                     y = .data[["Ref"]],
                     fill = .data[["Similarity"]])
    ) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_gradient2(
            low = "blue",
            high = "red",
            mid = "white",
            midpoint = 0,
            limit = color_limits,
            space = "Lab",
            name = paste(stringr::str_to_title(metric_name),
                         "\nSimilarity")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45,
                                                vjust = 1,
                                                size = 10,
                                                hjust = 1),
            axis.text.y = ggplot2::element_text(size = 10),
            plot.title = ggplot2::element_text(hjust = 0.5,
                                               size = 14),
            legend.title = ggplot2::element_text(size = 10),
            panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            x = "Query Dataset PCs",
            y = "Reference Dataset PCs",
            title = title_text
        )

    # Add similarity values if requested
    if (show_values) {
        if (!is.null(p_values) && show_significance) {
            # Use the pre-created Label column and add color mapping
            pc_plot <- pc_plot +
                ggplot2::geom_text(
                    ggplot2::aes(label = .data[["Label"]],
                                 color = .data[["Significant"]]),
                    size = 3, fontface = "bold"
                ) +
                ggplot2::scale_color_manual(
                    values = c("TRUE" = "black", "FALSE" = "gray60"),
                    guide = "none"
                )
        } else {
            # Use the simple Label column
            pc_plot <- pc_plot +
                ggplot2::geom_text(
                    ggplot2::aes(label = .data[["Label"]]),
                    size = 3, color = "black"
                )
        }
    }

    # Add significance legend if applicable
    if (!is.null(p_values) && show_significance && show_values) {
        pc_plot <- pc_plot +
            ggplot2::labs(
                caption = paste0("* indicates p < ",
                                 significance_threshold,
                                 " (based on ",
                                 x[["n_permutations"]],
                                 " permutations)")
            ) +
            ggplot2::theme(
                plot.caption = ggplot2::element_text(hjust = 0,
                                                     size = 8,
                                                     color = "gray50")
            )
    }

    return(pc_plot)
}
