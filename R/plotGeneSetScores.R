#' @title Visualization of gene sets or pathway scores on dimensional reduction plot
#' 
#' @description
#' Plot gene sets or pathway scores on PCA, TSNE, or UMAP. Single cells are color-coded by scores of gene sets or pathways.
#' 
#' @details 
#' This function plots gene set scores on reduced dimensions such as PCA, t-SNE, or UMAP. 
#' It extracts the reduced dimensions from the provided SingleCellExperiment object.
#' Gene set scores are visualized as a scatter plot with colors indicating the scores.
#' For PCA, the function automatically includes the percentage of variance explained 
#' in the plot's legend.
#'          
#' @param se_object An object of class \code{\linkS4class{SingleCellExperiment}} containing numeric expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param method A character string indicating the method for visualization ("PCA", "TSNE", or "UMAP").
#' @param score_col A character string representing the name of the score_col (score) in the colData(se_object) to plot.
#' @param pc_subset An optional vector specifying the principal components (PCs) to include in the plot if method = "PCA". 
#'        Default is 1:5.
#'
#' @return A ggplot2 object representing the gene set scores plotted on the specified reduced dimensions.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("query_data")
#' 
#' # Plot gene set scores on PCA
#' plotGeneSetScores(se_object = query_data, 
#'                   method = "PCA", 
#'                   score_col = "gene_set_scores",
#'                   pc_subset = 1:5)
#'                   
#' # Note: Users can provide their own gene set scores in the colData of the 'se_object' object, 
#' # using any dimension reduction of their choice.
#'
# Function to plot gene get scores in PCA space
plotGeneSetScores <- function(se_object, 
                              method = c("PCA", "TSNE", "UMAP"), 
                              score_col,
                              pc_subset = 1:5) {
    
    # Match method argument
    method <- match.arg(method)
    
    # Check if dimention reduction method is present in reference's reducedDims
    if (!(method %in% names(reducedDims(se_object)))) {
        stop("\'se_object\' must have pre-computed \'", method, "\' in \'reducedDims\'.")
    }
    
    # Check standard input arguments
    argumentCheck(query_data = se_object,
                  pc_subset_query = pc_subset)
    
    # Check if score_col is a valid column name in se_object
    if (!score_col %in% colnames(colData(se_object))) {
        stop("score_col: '", score_col, "' is not a valid column name in se_object.")
    }
    
    # Create the plot object
    if (method == "PCA") {
        
        # PCA data
        plot_mat <- reducedDim(se_object, "PCA")[, pc_subset]
        # Modify column names to include percentage of variance explained
        colnames(plot_mat) <- paste0(
            "PC", pc_subset, 
            " (", sprintf("%.1f%%", 
                          attributes(reducedDim(se_object, "PCA"))[["percentVar"]][pc_subset]), ")")
    } else if (method == "TSNE") {
        
        # TSNE data
        plot_mat <- reducedDim(se_object, "TSNE")
        
    } else if (method == "UMAP") {
        
        # UMAP data
        plot_mat <- reducedDim(se_object, "UMAP")
    }
    
    # Create all possible pairs of specified PCs
    plot_names <- colnames(plot_mat)
    pairs <- expand.grid(x = plot_names, y = plot_names)
    pairs <- pairs[pairs[["x"]] != pairs[["y"]], ]
    # Create a new data frame with all possible pairs of specified PCs
    data_pairs_list <- lapply(seq_len(nrow(pairs)), function(i) {
        x_col <- pairs[["x"]][i]
        y_col <- pairs[["y"]][i]
        data_frame <- data.frame(plot_mat[, c(x_col, y_col)])
        colnames(data_frame) <- c("x_value", "y_value")
        data_frame[["x"]] <- x_col
        data_frame[["y"]] <- y_col
        return(data_frame)
    })
    # Plot data
    data_pairs <- do.call(rbind, data_pairs_list)
    # Remove redundant data (to avoid duplicated plots)
    data_pairs <- data_pairs[as.numeric(data_pairs[["x"]]) < 
                                 as.numeric(data_pairs[["y"]]),]
    data_pairs[["Scores"]] <- se_object[[score_col]]
    # Create the ggplot object (with facets if PCA)
    plot_obj <- ggplot2::ggplot(
        data_pairs, ggplot2::aes(x = .data[["x_value"]], 
                                 y = .data[["y_value"]], 
                                 color = .data[["Scores"]])) +
        ggplot2::geom_point(size = 1, alpha = 0.5) + 
        ggplot2::xlab("") + ggplot2::ylab("") + 
        ggplot2::scale_color_gradientn(
            colors = c("#2171B5", "#8AABC1", "#FFEDA0", "#E6550D"), 
            values = seq(0, 1, by = 1/3), 
            limits = c(0, max(data_pairs$Scores))) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(.data[["y"]]), 
            cols = ggplot2::vars(.data[["x"]]), scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray", 
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14, face = "bold", 
                                               hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12), 
            axis.text = ggplot2::element_text(size = 10))
    return(plot_obj)
}
