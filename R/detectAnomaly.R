#' 
#' @importFrom methods is
#' @importFrom stats na.omit predict qnorm
#' @importFrom utils tail
#' 
#' @title PCA Anomaly Scores via Isolation Forests with Visualization
#'
#' @description 
#' This function detects anomalies in single-cell data by projecting the data onto a PCA space and using an isolation forest 
#' algorithm to identify anomalies.
#'
#' @details This function projects the query data onto the PCA space of the reference data. An isolation forest is then built on the 
#' reference data to identify anomalies in the query data based on their PCA projections. If no query dataset is provided by the user,
#' the anomaly scores are computed on the reference data itself. Anomaly scores for the data with all combined cell types are also
#' provided as part of the output.
#' 
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data An optional \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells. 
#' If NULL, then the isolation forest anomaly scores are computed for the reference data. Default is NULL.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param n_components An integer specifying the number of principal components to use. Default is 10.
#' @param n_tree An integer specifying the number of trees for the isolation forest. Default is 500
#' @param anomaly_treshold A numeric value specifying the threshold for identifying anomalies, Default is 0.5.
#' @param ... Additional arguments passed to the `isolation.forest` function.
#' 
#' @return A list containing the following components for each cell type and the combined data:
#' \item{anomaly_scores}{Anomaly scores for each cell in the query data.}
#' \item{anomaly}{Logical vector indicating whether each cell is classified as an anomaly.}
#' \item{reference_mat_subset}{PCA projections of the reference data.}
#' \item{query_mat_subset}{PCA projections of the query data (if provided).}
#' \item{var_explained}{Proportion of variance explained by the retained principal components.}
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
#' plot(anomaly_output, cell_type = "CD8", pc_subset = c(1:5), data_type = "query")
#' 
# Function to perform diagnostics using isolation forest with PCA and visualization
detectAnomaly <- function(reference_data, 
                          query_data = NULL, 
                          ref_cell_type_col,
                          query_cell_type_col, 
                          n_components = 10,
                          n_tree = 500,
                          anomaly_treshold = 0.5,
                          ...) {
    
    # Check whether the anlaysis is done only for one dataset
  if (is.null(query_data)) {
      include_query_in_output <- FALSE
  } else{
      if(is.null(query_cell_type_col))
          stop("If \'query_data\' is not NULL, a value for \'query_cell_type_col\' must be provided.")
      include_query_in_output <- TRUE
  }
    
  if(!is.null(n_components)){
      reference_mat <- reducedDim(reference_data, "PCA")[, 1:n_components]
      if(include_query_in_output){
          # Get PCA data from reference and query datasets (query data projected onto PCA space of reference dataset)
          pca_output <- projectPCA(query_data = query_data, reference_data = reference_data, 
                                   query_cell_type_col = query_cell_type_col, ref_cell_type_col = ref_cell_type_col,
                                   n_components = n_components, return_value = "list")
          query_mat <- pca_output$query[, paste0("PC", 1:n_components)]
      }
  } else{
      reference_mat <- t(as.matrix(assay(reference_data, "logcounts")))
      if(include_query_in_output){
          query_mat <- t(as.matrix(assay(query_data, "logcounts")))
      }
  }
  
  # List to store output
  output <- list()
  
  # Extract reference and query annotations
  reference_labels <- reference_data[[ref_cell_type_col]]
  if(!include_query_in_output){
      cell_types <- c(as.list(na.omit(unique(reference_labels))),
                      list(na.omit(unique(reference_labels))))
  } else{
      query_labels <- query_data[[query_cell_type_col]]
      cell_types <- c(as.list(na.omit(intersect(unique(reference_labels), unique(query_labels)))),
                      list(na.omit(intersect(unique(reference_labels), unique(query_labels)))))
  }

  for (cell_type in cell_types) {
    
    # Filter reference and query PCA data for the current cell type
    reference_mat_subset <- na.omit(reference_mat[reference_labels %in% cell_type,])
    
    # Build isolation forest on reference PCA data for this cell type
    isolation_forest <- isotree::isolation.forest(reference_mat_subset, ntree = n_tree, ...)
      
    # Calculate anomaly scores for query data (scaled by reference path length)
    reference_anomaly_scores <- predict(isolation_forest, newdata = reference_mat_subset, type = "score")
    if(include_query_in_output){
        query_mat_subset <- na.omit(query_mat[query_labels %in% cell_type,])
        query_anomaly_scores <- predict(isolation_forest, newdata = query_mat_subset, type = "score")
    }

    # Store cell type anomaly scores and PCA data
    list_name <- ifelse(length(cell_type) == 1, cell_type, "Combined")
    output[[list_name]] <- list()
    output[[list_name]]$reference_anomaly_scores <- reference_anomaly_scores
    output[[list_name]]$reference_anomaly <- reference_anomaly_scores > anomaly_treshold
    output[[list_name]]$reference_mat_subset <- reference_mat_subset
    if(include_query_in_output){
        output[[list_name]]$query_mat_subset <- query_mat_subset
        output[[list_name]]$query_anomaly_scores <- query_anomaly_scores
        output[[list_name]]$query_anomaly <- query_anomaly_scores > anomaly_treshold
    }
    if(!is.null(n_components))
        output[[list_name]]$var_explained <- (attributes(reducedDim(reference_data, "PCA"))$varExplained[1:n_components]) /
        sum(attributes(reducedDim(reference_data, "PCA"))$varExplained) 
  }
  
  # Set the class of the output
  class(output) <- c(class(output), "detectAnomaly")
  
  # Return anomaly, PCA data and optional PCA anomaly plots for each cell type
  return(output)
}
