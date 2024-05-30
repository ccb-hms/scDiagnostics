#' 
#' @importFrom methods is
#' @importFrom stats na.omit predict qnorm
#' @importFrom utils tail
#' 
#' @title PCA Anomaly Scores via Isolation Forests with Visualization
#'
#' @description \code{detectAnomaly} performs diagnostics using isolation forest with PCA and visualization. 
#' It takes reference and query \code{\linkS4class{SingleCellExperiment}} objects, their corresponding labels, and various parameters to perform 
#' the analysis. The function returns a list containing the results for each cell type, including anomaly scores, anomaly IDs, 
#' PCA data, and optional PCA anomaly plots.
#'
#' @details
#' This function first applies PCA to the entire reference expression data and predicts PCA scores for the query data if it is provided. 
#' It builds isolation forests and performs diagnostics for each cell type, calculating anomaly scores for the query data (or reference 
#' data if the query data is not provided).
#' Isolation Forest is an algorithm for anomaly detection that works by building an ensemble of isolation trees. It is based on the 
#' idea that anomalies are more susceptible to isolation than normal instances.
#' The part where we project the query data onto the PCA space of the reference data is done by using the `predict` function on the PCA model with the query expression data. This allows us to transform the query data into the same PCA space as the reference data, which is necessary for the isolation forest analysis.
#' 
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data An optional \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells. 
#' If NULL, then the isolation forest anomaly scores are computed for the reference data. Default is NULL.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param n_components An integer specifying the number of principal components to use. Default is 10.
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 500
#' @param anomaly_treshold A numeric value specifying the threshold for identifying anomalies, Default is 0.5.
#' @param ... Additional arguments passed to the `isolation.forest` function.
#' 
#' @return A list containing the results for each cell type, including anomaly scores, anomaly IDs and PCA data.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.detectAnomaly}}
#' 
#' @examples
#' # Load required libraries
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scran)
#' library(scater)
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
#' # log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#' 
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' 
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- getTopHVGs(ref_data, n = 2000)
#' query_var <- getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#' 
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset)
#' 
#' # Store PCA anomaly data and plots
#' anomaly_output <- detectAnomaly(ref_data_subset, query_data_subset,
#'                                 ref_cell_type_col = "reclustered.broad", 
#'                                 query_cell_type_col = "labels",
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5) 
#' 
#' # Plot the output for a cell type
#' plot(anomaly_output, cell_type = "CD8", pc_subset = c(1:5))
#' 
#' 
# Function to perform diagnostics using isolation forest with PCA and visualization
detectAnomaly <- function(reference_data, 
                          query_data = NULL, 
                          query_cell_type_col, 
                          ref_cell_type_col, 
                          n_components = 10,
                          n_tree = 500,
                          anomaly_treshold = 0.5,
                          ...) {
    
    # Check whether the anlaysis is done only for one dataset
  if (is.null(query_data)) {
      query_data <- reference_data
      query_cell_type_col <- ref_cell_type_col
  }
    
  if(!is.null(n_components)){
      # Get PCA data from reference and query datasets (query data projected onto PCA space of reference dataset)
      pca_output <- projectPCA(query_data = query_data, reference_data = reference_data, 
                               n_components = n_components, return_value = "list")
      reference_mat <- pca_output$ref[, paste0("PC", 1:n_components)]
      query_mat <- pca_output$query[, paste0("PC", 1:n_components)]
  } else{
      reference_mat <- t(as.matrix(assay(reference_data, "logcounts")))
      query_mat <- t(as.matrix(assay(query_data, "logcounts")))
  }
  
  # List to store output
  output <- list()
  
  # Extract reference and query annotations
  reference_labels <- reference_data[[ref_cell_type_col]]
  query_labels <- query_data[[query_cell_type_col]]
  
  # Build isolation forests and perform diagnostics for each cell type
  cell_types <- na.omit(intersect(unique(reference_labels), unique(query_labels)))
  for (cell_type in cell_types) {
    
    # Filter reference and query PCA data for the current cell type
    reference_mat_subset <- na.omit(reference_mat[reference_labels == cell_type,])
    query_mat_subset <- na.omit(query_mat[query_labels == cell_type,])
    
    # Build isolation forest on reference PCA data for this cell type
    isolation_forest <- isotree::isolation.forest(reference_mat_subset, ntree = n_tree, ...)
      
    # Calculate anomaly scores for query data (scaled by reference path length)
    anomaly_scores <- predict(isolation_forest, newdata = query_mat_subset, type = "score")

    # Create list of output for cell type
    output[[paste0(cell_type)]] <- list()
    
    # Store cell type anomaly scores and PCA data
    output[[paste0(cell_type)]]$anomaly_scores <- anomaly_scores
    output[[paste0(cell_type)]]$anomaly <- anomaly_scores > anomaly_treshold
    output[[paste0(cell_type)]]$reference_mat_subset <- reference_mat_subset
    output[[paste0(cell_type)]]$query_mat_subset <- query_mat_subset
    if(!is.null(n_components))
        output[[paste0(cell_type)]]$var_explained <- (attributes(reducedDim(reference_data, "PCA"))$varExplained[1:n_components]) /
        sum(attributes(reducedDim(reference_data, "PCA"))$varExplained) 
  }
  
  # Set the class of the output
  class(output) <- c(class(output), "detectAnomaly")
  
  # Return anomaly, PCA data and optional PCA anomaly plots for each cell type
  return(output)
}
