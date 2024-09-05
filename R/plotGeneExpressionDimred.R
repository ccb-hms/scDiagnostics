#' @title Visualize gene expression on a dimensional reduction plot
#'
#' @description
#' This function plots gene expression on a dimensional reduction plot using methods like t-SNE, UMAP, or PCA. Each single cell is color-coded based on the expression of a specific gene or feature.
#'
#' @param se_object An object of class \code{\linkS4class{SingleCellExperiment}} containing log-transformed expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param method The reduction method to use for visualization. It should be one of the supported methods: "TSNE", "UMAP", or "PCA".
#' @param pc_subset An optional vector specifying the principal components (PCs) to include in the plot if method = "PCA". 
#'        Default is 1:5.
#' @param feature A character string representing the name of the gene or feature to be visualized.
#'
#' @importFrom SummarizedExperiment assay
#' @import SingleCellExperiment
#'
#' @return A ggplot object representing the dimensional reduction plot with gene expression.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("query_data")
#' 
#' # Plot gene expression on PCA plot
#' plotGeneExpressionDimred(se_object = query_data, 
#'                          method = "PCA", 
#'                          pc_subset = 1:5, 
#'                          feature = "VPREB3")
#' 
# Function to plot gene expression in the reduced dimension
plotGeneExpressionDimred <- function(se_object, 
                                     method = c("TSNE", "UMAP", "PCA"), 
                                     pc_subset = 1:5, 
                                     feature) {
    
    # Check standard input arguments
    argumentCheck(query_data = se_object,
                  pc_subset_query = pc_subset)
    
    # Match argument for method
    method <- match.arg(method)
    
    # Check if feature is available
    if (!feature %in% rownames(assay(se_object, "logcounts"))) {
        stop("Specified feature does not exist in the expression matrix.")
    }
    
    # Extract gene expression vector
    expression <- assay(se_object, "logcounts")[feature, ]
    
    if(method %in% c("TSNE", "UMAP")){
        
        # Extract dimension reduction coordinates from SingleCellExperiment object
        reduction <- reducedDim(se_object, method)
        
        # Prepare data for plotting
        df <- data.frame(Dim1 = reduction[, 1], Dim2 = reduction[, 2], 
                         Expression = expression)
        
        # Create the plot object
        plot_obj <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Dim1"]], y = .data[["Dim2"]])) +
            ggplot2::geom_point(ggplot2::aes(color = .data[["Expression"]])) +
            ggplot2::scale_color_gradient(low = "grey90", high = "blue") +
            ggplot2::xlab("Dimension 1") +
            ggplot2::ylab("Dimension 2") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                strip.background = ggplot2::element_rect(fill = "grey85", 
                                                         color = "grey70"),   
                strip.text = ggplot2::element_text(size = 10, 
                                                   face = "bold", 
                                                   color = "black"), 
                axis.title = ggplot2::element_blank(), 
                axis.text = ggplot2::element_text(size = 10), 
                panel.grid = ggplot2::element_blank(), 
                panel.background = ggplot2::element_rect(fill = "white", 
                                                         color = "black"), 
                legend.position = "right", 
                plot.title = ggplot2::element_text(size = 14, hjust = 0.5), 
                plot.background = ggplot2::element_rect(fill = "white")) 
        
    } else if (method == "PCA"){
        
        # Check input for pc_subset
        if(!all(pc_subset %in% seq_len(ncol(reducedDim(se_object, "PCA")))))
            stop("\'pc_subset\' is out of range.")
        
        # PCA data
        plot_mat <- reducedDim(se_object, "PCA")[, pc_subset]
        # Modify column names to include percentage of variance explained
        colnames(plot_mat) <- paste0("PC", pc_subset, 
                                     " (", sprintf("%.1f%%", 
                                                   attributes(reducedDim(se_object, "PCA"))[["percentVar"]][pc_subset]), ")")
        
        # Create all possible pairs of specified PCs
        plot_names <- colnames(plot_mat)
        pairs <- expand.grid(x = plot_names, y = plot_names)
        pairs <- pairs[pairs$x != pairs$y, ]
        # Create a new data frame with all possible pairs of specified PCs
        data_pairs_list <- lapply(seq_len(nrow(pairs)), function(i) {
            x_col <- pairs$x[i]
            y_col <- pairs$y[i]
            data_frame <- data.frame(plot_mat[, c(x_col, y_col)])
            colnames(data_frame) <- c("x_value", "y_value")
            data_frame$x <- x_col
            data_frame$y <- y_col
            data_frame
        })
        # Plot data
        data_pairs <- do.call(rbind, data_pairs_list)
        # Remove redundant data (to avoid duplicated plots)
        data_pairs <- data_pairs[as.numeric(data_pairs$x) < 
                                     as.numeric(data_pairs$y),]
        data_pairs$Expression <- expression
        
        # Create the ggplot object (with facets if PCA)
        plot_obj <- ggplot2::ggplot(
            data_pairs, ggplot2::aes(x = .data[["x_value"]], 
                                     y = .data[["y_value"]], 
                                     color = .data[["Expression"]])) +
            ggplot2::geom_point(size = 1, alpha = 0.5) + 
            ggplot2::xlab("") + ggplot2::ylab("") + 
            ggplot2::scale_color_gradient(low = "grey90", high = "blue") +
            ggplot2::facet_grid(rows = ggplot2::vars(.data[["y"]]), 
                                cols = ggplot2::vars(.data[["x"]]), 
                                scales = "free") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_line(color = "gray", 
                                                         linetype = "dotted"),
                plot.title = ggplot2::element_text(size = 14, face = "bold", 
                                                   hjust = 0.5),
                axis.title = ggplot2::element_text(size = 12), 
                axis.text = ggplot2::element_text(size = 10))
    }
    return(plot_obj)
}
