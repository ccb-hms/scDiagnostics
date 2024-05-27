#' @title Compute Distance Diagnostics Between Reference and Query Data
#'
#' @description This function computes the distances within the reference dataset and the distances from each query sample to all 
#' reference samples for each cell type. It uses PCA for dimensionality reduction and Euclidean distance for distance calculation.
#'
#' @details The function first performs PCA on the reference dataset and projects the query dataset onto the same PCA space. 
#' It then computes pairwise Euclidean distances within the reference dataset for each cell type, as well as distances from each 
#' query sample to all reference samples of the same cell type. The results are stored in a list, with one entry per cell type.
#'
#' @param query_data A SingleCellExperiment object containing the data to be projected.
#' @param reference_data A SingleCellExperiment object containing the reference data with pre-computed PCA.
#' @param n_components An integer specifying the number of principal components to use for projection. Defaults to 10. 
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#'
#' @return A list containing distance data for each cell type. Each entry in the list contains:
#' \describe{
#'   \item{ref_distances}{A vector of all pairwise distances within the reference subset for the cell type.}
#'   \item{query_to_ref_distances}{A matrix of distances from each query sample to all reference samples for the cell type.}
#' }
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load required libraries
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scran)
#' library(scater)
#'
#' # Load data (replace with your data loading)
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
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
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR::SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' 
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- scran::getTopHVGs(ref_data, n = 2000)
#' query_var <- scran::getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#' 
#' # Plot the PC data
#' distance_data <- computeDistanceDiagnostics(query_data, reference_data, 
#'                                             n_components = 10, 
#'                                             query_cell_type_col = "labels", 
#'                                             ref_cell_type_col = "reclustered.broad",
#'                                             pc_subset = c(1:10)) 
#' 
#' # Identify outliers for CD4
#' cd4_anomalites <- detectAnomaly(reference_data = ref_data_subset, query_data = query_data_subset, 
#'                                 query_cell_type_col = "labels", 
#'                                 ref_cell_type_col = "reclustered.broad",
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5,
#'                                 store_plots = TRUE)$CD4
#' cd4_top_anomaly <- names(which.max(cd4_anomalites$anomaly_scores))
#' 
#' # Plot the densities of the distances
#' plot(distance_data, cell_type = "CD4", sample_name = cd4_top_anomaly)
#' 
# Function to compute distances from 
computeDistanceDiagnostics <- function(query_data, reference_data, 
                                       query_cell_type_col, 
                                       ref_cell_type_col,
                                       n_components = 10, 
                                       pc_subset = c(1:5)) {
    
    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data, reference_data = reference_data, 
                             n_components = n_components, 
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col, 
                             return_value = "list")
    
    # Get unique cell types
    unique_cell_types <- na.omit(unique(c(colData(reference_data)[[ref_cell_type_col]],
                                          colData(query_data)[[query_cell_type_col]])))
    
    # Create a list to store distance data for each cell type
    distance_data <- list()
    
    # Function to compute Euclidean distance between a vector and each row of a matrix
    .compute_distances <- function(matrix, vector) {

        # Apply the distance function to each row of the matrix
        distances <- apply(matrix, 1, function(row) {
            sqrt(sum((row - vector) ^ 2))
        })
        
        return(distances)
    }
    
    for (cell_type in unique_cell_types) {
        
        # Subset principal component scores for current cell type
        ref_subset_scores <- pca_output$ref[which(cell_type == reference_data[[ref_cell_type_col]]), pc_subset]
        query_subset_scores <- pca_output$query[which(cell_type == query_data[[query_cell_type_col]]), pc_subset]
        
        # Compute all pairwise distances within the reference subset
        ref_distances <- as.vector(dist(ref_subset_scores))
        
        # Compute distances from each query sample to all reference samples
        query_to_ref_distances <- apply(query_subset_scores, 1, function(query_sample, ref_subset_scores) {
            .compute_distances(ref_subset_scores, query_sample)
        }, ref_subset_scores = ref_subset_scores)
        
        # Store the distances
        distance_data[[cell_type]] <- list(
            ref_distances = ref_distances,
            query_to_ref_distances = t(query_to_ref_distances)
        )
    }
    
    # Add class of object
    class(distance_data) <- c(class(distance_data), "computeDistanceDiagnostics")
    
    # Return the distance data
    return(distance_data)
}