#' 
#' @title PCA Anomaly Scores via Isolation Forests with Visualization
#'
#' @description \code{diagnose_iForest} performs diagnostics using isolation forest with PCA and visualization. 
#' It takes reference and query SingleCellExperiment objects, their corresponding labels, and various parameters to perform the analysis. 
#' The function returns a list containing the results for each cell type, including anomaly scores, outlier IDs, PCA data, 
#' and optional PCA anomaly plots.
#'
#' @details
#' This function first gets common genes between the reference and query data, and then applies log normalization to the expression data. 
#' It selects highly variable genes from the reference data and extracts the expression data and labels.
#' The function then applies PCA to the entire reference expression data and predicts PCA scores for the query data. 
#' It builds isolation forests and performs diagnostics for each cell type, calculating anomaly scores for the query data.
#' The function also provides optional analysis of anomaly scores and creates a list of output for each cell type, 
#' including PCA anomaly plots if `store_plots` is set to TRUE. 
#' Finally, it returns a list containing the results for each cell type.
#' Isolation Forest is an algorithm for anomaly detection that works by building an ensemble of isolation trees. It is based on the 
#' idea that anomalies are more susceptible to isolation than normal instances.
#' The part where we project the query data onto the PCA space of the reference data is done by using the `predict` function on the PCA model with the query expression data. This allows us to transform the query data into the same PCA space as the reference data, which is necessary for the isolation forest analysis.
#' 
#' @param reference_sce A SingleCellExperiment object containing the reference data.
#' @param query_sce A SingleCellExperiment object containing the query data.
#' @param reference_labels A vector of labels for the reference data.
#' @param query_labels A vector of labels for the query data.
#' @param n_pca An integer specifying the number of principal components to use. Default is 10.
#' @param n_hvg An integer specifying the number of highly variable genes to use. Default is 2000.
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 100.
#' @param anomaly_treshold A numeric value specifying the threshold for identifying outliers. Default is 0.5.
#' @param verbose A logical value indicating whether to print verbose output. Default is FALSE.
#' @param store_plots A logical value indicating whether to store PCA anomaly plots. Default is TRUE.
#' @param ... Additional arguments passed to the `isolation.forest` function.
#' 
#' @return A list containing the results for each cell type, including anomaly scores, outlier IDs, PCA data, and 
#' optional PCA anomaly plots.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @examples
#' # Load data
#' sce <- scRNAseq::HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # log transform datasets
#' ref_data <- scuttle::logNormCounts(ref_data)
#' query_data <- scuttle::logNormCounts(query_data)
#'
#' # Cell type annotation via SingleR (or any other cell type annotation method)
#' scores <- SingleR::SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#'
#' # Store PCA anomaly data and plots
#' iForest_output <- diagnose_iForest(reference_sce = ref_data, query_sce = query_data, 
#'                                    reference_labels = ref_data$reclustered.broad, query_labels = query_data$labels,
#'                                    n_pca = 10,
#'                                    n_hvg = 2000,
#'                                    n_tree = 100,
#'                                    anomaly_treshold = 0.5,
#'                                    verbose = FALSE,
#'                                    store_plots = TRUE)
#'
#' # Plot results for CD4
#' plot_list <- list(iForest_output$CD4$PC_plots$PC1_PC2, iForest_output$CD4$PC_plots$PC3_PC4, 
#'                   iForest_output$CD4$PC_plots$PC5_PC6, iForest_output$CD4$PC_plots$PC7_PC8)
#' gridExtra::grid.arrange(grobs = plot_list, ncol = 2)#' 
#' 
# Function to perform diagnostics using isolation forest with PCA and visualization
diagnose_iForest <- function(reference_sce, query_sce, 
                             reference_labels, query_labels, 
                             n_pca = 10,
                             n_hvg = 2000,
                             n_tree = 100,
                             anomaly_treshold = 0.5,
                             verbose = FALSE,
                             store_plots = TRUE,
                             ...) {
  
  # Get common genes
  common_genes <- BiocGenerics::intersect(rownames(reference_sce), rownames(query_sce))
  # Keep common genes
  reference_sce <- reference_sce[common_genes,]
  query_sce <- query_sce[common_genes,]
  
  # Get log counts
  reference_sce <- scuttle::logNormCounts(reference_sce)
  query_sce <- scuttle::logNormCounts(query_sce)
  
  # Selecting highly variable genes (from reference data)
  ref_hvg <- scran::getTopHVGs(reference_sce, n = n_hvg)
  
  # Extract expression data and labels (assuming cell type is included)
  reference_expr <- t(as.matrix(SummarizedExperiment::assay(reference_sce[ref_hvg,], "logcounts")))
  query_expr <- t(as.matrix(SummarizedExperiment::assay(query_sce[ref_hvg,], "logcounts")))

  # Apply PCA to the entire reference expression data
  pca_model <- stats::prcomp(reference_expr, scale. = TRUE, center = TRUE, rank. = n_pca)
  reference_pca <- pca_model$x[, 1:n_pca]
  
  # Predict PCA scores for the query data
  query_pca <- predict(pca_model, query_expr)[, 1:n_pca]
  
  # List to store output
  output <- list()
  
  # Build isolation forests and perform diagnostics for each cell type
  cell_types <- BiocGenerics::unique(query_labels)[!is.na(BiocGenerics::unique(query_labels))]  # Extract unique cell types from query data
  for (cell_type in cell_types) {
    
    # Filter reference and query PCA data for the current cell type
    reference_pca_subset <- na.omit(reference_pca[reference_labels == cell_type,])
    query_pca_subset <- na.omit(query_pca[query_labels == cell_type,])
    
    # Build isolation forest on reference PCA data for this cell type
    isolation_forest <- isotree::isolation.forest(reference_pca_subset, ntree = n_tree, ...)
      
    # Calculate anomaly scores for query data (scaled by reference path length)
    anomaly_scores <- predict(isolation_forest, newdata = query_pca_subset, type = "score")

    # Optional analysis of anomaly scores
    if(verbose){
      outlier_proportion <- sum(anomaly_scores > anomaly_treshold) / nrow(query_pca_subset)
      cat("Proportion of potential outliers:", outlier_proportion, "\n")
    }

    # Create list of output for cell type
    output[[paste0(cell_type)]] <- list()
    
    # Visualization of anomaly scores for PCs
    if(store_plots){
      
      output[[paste0(cell_type)]][["PC_plots"]] <- list()
      plot_combn <-  t(utils::combn(n_pca, 2))
      
      # Looping over all combinations of PCs
      for(plot_id in 1:nrow(plot_combn)){
        
        output[[paste0(cell_type)]][["PC_plots"]][[paste0("PC", plot_combn[plot_id, 1], "_PC", plot_combn[plot_id, 2])]] <- 
          ggplot2::ggplot(data.frame(x = query_pca_subset[, plot_combn[plot_id, 1]],
                                     y = query_pca_subset[, plot_combn[plot_id, 2]],
                                     anomaly_score = anomaly_scores),
                          ggplot2::aes(x = x, y = y, color = anomaly_score)) +
          ggplot2::geom_point(size = 4, alpha = 0.7, shape = 19) +
          ggplot2::labs(title = paste0("PCA colored by Anomaly Scores (Isolation Forest - ", cell_type, ")"),
                        x = paste0("PC", plot_combn[plot_id, 1]), y = paste0("PC", plot_combn[plot_id, 2]), color = "Anomaly Score") +
          ggplot2::scale_color_gradient(low = "green", high = "red", name = "Anomaly Score") +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, color = "black"),
                         legend.position = "bottom",
                         legend.title = ggplot2::element_text(face = "bold", size = 14, color = "black"),
                         legend.text = ggplot2::element_text(size = 12, color = "black"),
                         panel.grid.major = ggplot2::element_line(color = "grey90", size = 0.2),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.border = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(color = "black"),
                         axis.title = ggplot2::element_text(face = "bold", size = 14, color = "black"),
                         axis.text = ggplot2::element_text(size = 12, color = "black"),
                         plot.background = ggplot2::element_rect(fill = "white", color = NA),
                         panel.background = ggplot2::element_rect(fill = "white", color = NA),
                         plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"))
      }
    }
    
    # Store cell type anomaly scores and PCA data
    output[[paste0(cell_type)]]$anomaly_scores <- anomaly_scores
    output[[paste0(cell_type)]]$outlier_id <- rownames(query_pca_subset)[which(anomaly_scores > anomaly_treshold)]
    output[[paste0(cell_type)]]$reference_pca_subset <- reference_pca_subset
    output[[paste0(cell_type)]]$query_pca_subset <- query_pca_subset
  }
  
  # Return anomaly, PCA data and optional PCA anomaly plots for each cell type
  return(output)
}
