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
#' @param se_object An object of class "SingleCellExperiment" containing numeric expression matrix and other metadata.
#'        It can be either a reference or query dataset.
#' @param method A character string indicating the method for visualization ("PCA", "TSNE", or "UMAP").
#' @param feature A character string representing the name of the feature (score) in the colData(query_data) to plot.
#' @param pc_subset An optional vector specifying the principal components (PCs) to include in the plot if method = "PCA". 
#'        Default is c(1:5).
#'
#' @return A ggplot2 object representing the gene set scores plotted on the specified reduced dimensions.
#' @export
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(AUCell)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' ## log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Run PCA on the query data
#' query_data <- runPCA(query_data)
#' 
#' # Compute scores using AUCell
#' expression_matrix <- assay(query_data, "logcounts")
#' cells_rankings <- AUCell_buildRankings(expression_matrix, plotStats = FALSE)
#' # Generate gene sets
#' gene_set1 <- sample(rownames(expression_matrix), 10)
#' gene_set2 <- sample(rownames(expression_matrix), 20)
#' gene_sets <- list(geneSet1 = gene_set1, geneSet2 = gene_set2)
#' 
#' # Calculate AUC scores for gene sets
#' cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)
#' 
#' # Assign scores to colData (users should ensure that the scores are present in the colData)
#' colData(query_data)$geneSetScores <- assay(cells_AUC)["geneSet1", ] 
#'
#' # Plot gene set scores on PCA
#' plotGeneSetScores(se_object = query_data, 
#'                   method = "PCA", 
#'                   feature = "geneSetScores",
#'                   pc_subset = c(1:5))
#'
#' # Note: Users can provide their own gene set scores in the colData of the 'se_object' object, 
#' # using any method of their choice.
#'
plotGeneSetScores <- function(se_object, 
                              method, 
                              feature,
                              pc_subset = c(1:5)) {

  # Check if the specified method is valid
  valid_methods <- c("PCA", "TSNE", "UMAP")
  if (!(method %in% valid_methods)) {
    stop("Invalid method. Please choose one of: ", paste(valid_methods, collapse = ", "))
  }

  # Create the plot object
  if (method == "PCA") {
      # Check if "PCA" is present in reference's reduced dimensions
      if (!"PCA" %in% names(reducedDims(se_object))) {
          stop("Reference data must have pre-computed PCA in \'reducedDims\'.")
      }
      
      # Check input for pc_subset
      if(!all(pc_subset %in% 1:ncol(reducedDim(se_object, "PCA"))))
          stop("\'pc_subset\' is out of range.")
      
      # PCA data
      plot_mat <- reducedDim(se_object, "PCA")[, pc_subset]
      # Modify column names to include percentage of variance explained
      colnames(plot_mat) <- paste0("PC", pc_subset, 
                                      " (", sprintf("%.1f%%", attributes(reducedDim(se_object, "PCA"))$varExplained[pc_subset] /
                                                        sum(attributes(reducedDim(se_object, "PCA"))$varExplained) * 100), ")")
  } else if (method == "TSNE") {
      # Check if "TSNE" is present in reference's reduced dimensions
      if (!"TSNE" %in% names(reducedDims(se_object))) {
          stop("Reference data must have pre-computed t-SNE in \'reducedDims\'.")
      }
      # TSNE data
      plot_mat <- reducedDim(se_object, "TSNE")
  } else if (method == "UMAP") {
      # Check if "UMAP" is present in reference's reduced dimensions
      if (!"UMAP" %in% names(reducedDims(se_object))) {
          stop("Reference data must have pre-computed UMAP in \'reducedDims\'.")
      }
      # UMAP data
      plot_mat <- reducedDim(se_object, "UMAP")
  }
  
  # Create all possible pairs of specified PCs
  plot_names <- colnames(plot_mat)
  pairs <- expand.grid(x = plot_names, y = plot_names)
  pairs <- pairs[pairs$x != pairs$y, ]
  # Create a new data frame with all possible pairs of specified PCs
  data_pairs_list <- lapply(1:nrow(pairs), function(i) {
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
  data_pairs <- data_pairs[as.numeric(data_pairs$x) < as.numeric(data_pairs$y),]
  data_pairs$Scores <- se_object[["geneSetScores"]]
  # Create the ggplot object (with facets if PCA)
  plot_obj <- ggplot2::ggplot(data_pairs, ggplot2::aes(x = x_value, y = y_value, color = Scores)) +
      ggplot2::geom_point(size = 1, alpha = 0.5) + 
      ggplot2::xlab("") + ggplot2::ylab("") + 
      ggplot2::scale_color_gradientn(colors = c("#2171B5", "#8AABC1", "#FFEDA0", "#E6550D"), 
                                     values = seq(0, 1, by = 1/3), 
                                     limits = c(0, max(data_pairs$Scores))) +
      ggplot2::facet_grid(rows = ggplot2::vars(y), cols = ggplot2::vars(x), scales = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                     plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10))
  return(plot_obj)
}
