#' @title Function to compute Bhattacharyya coefficients and Hellinger distances
#'
#' @description This function computes Bhattacharyya coefficients and Hellinger distances to quantify the overlap of density 
#' distributions between query samples and reference data for each cell type.

#'
#' @details This function first computes distance data using the \code{calculateDistanceDiagnostics} function, which calculates 
#' pairwise distances between samples within the reference data and between query samples and reference samples in the PCA space.
#' Bhattacharyya coefficients and Hellinger distances are calculated to quantify the overlap of density distributions between query 
#' samples and reference data for each cell type. Bhattacharyya coefficient measures the similarity of two probability distributions, 
#' while Hellinger distance measures the distance between two probability distributions.
#'
#' Bhattacharyya coefficients range between 0 and 1. A value closer to 1 indicates higher overlap between distributions, while a value 
#' closer to 0 indicates lower overlap.
#'
#' Hellinger distances range between 0 and 1. A value closer to 0 indicates higher similarity between distributions, while a value 
#' closer to 1 indicates lower similarity.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param sample_names A character vector specifying the names of the query samples for which to compute distance measures.
#' @param n_components An integer specifying the number of principal components to use for projection. Defaults to 10. 
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
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset, ncomponents = 50)
#' 
#' # Plot the PC data
#' distance_data <- calculateDistanceDiagnostics(query_data, reference_data, 
#'                                               n_components = 10, 
#'                                               query_cell_type_col = "labels", 
#'                                               ref_cell_type_col = "reclustered.broad",
#'                                               pc_subset = c(1:10)) 
#' 
#' # Identify outliers for CD4
#' cd4_anomalites <- detectAnomaly(query_data_subset, ref_data_subset, 
#'                                 query_cell_type_col = "labels", 
#'                                 ref_cell_type_col = "reclustered.broad",
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5,
#'                                 store_plots = TRUE)$CD4
#' cd4_top_anomaly <- names(which.max(cd4_anomalites$anomaly_scores))
#' 
#' # Get overlap measures
#' overlap_measures <- distanceDensityOverlapMeasures(query_data_subset,ref_data_subset, 
#'                                                    sample_names = cd4_top5_anomalies,
#'                                                    n_components = 10, 
#'                                                    query_cell_type_col = "labels", 
#'                                                    ref_cell_type_col = "reclustered.broad",
#'                                                    pc_subset = c(1:10))
#' 
# Function to compute Bhattacharyya coefficients and Hellinger distances
distanceDensityOverlapMeasures <- function(query_data, reference_data, 
                                           query_cell_type_col, 
                                           ref_cell_type_col,
                                           sample_names,
                                           n_components = 10, 
                                           pc_subset = c(1:5)) {

    # Check if samples are available in data for that cell type
    if(!all(sample_names %in% colnames(query_data)))
        stop("One or more specified 'sample_names' are not available for that cell type.")
    
    # Compute distance data
    query_data_subset <- query_data[, sample_names]
    distance_data <- calculateDistanceDiagnostics(query_data = query_data_subset, reference_data = reference_data, 
                                                  query_cell_type_col = query_cell_type_col, 
                                                  ref_cell_type_col = ref_cell_type_col,
                                                  n_components = n_components, 
                                                  pc_subset = pc_subset)
    
    # Initialize empty lists to store results
    bhattacharyya_list <- list()
    hellinger_list <- list()
    
    # Iterate over each cell type
    for (cell_type in names(distance_data)) {
        
        # Extract distances within the reference dataset for the current cell type
        ref_distances <- distance_data[[cell_type]]$ref_distances
        
        # Compute density of reference distances
        ref_density <- density(ref_distances)
        
        # Initialize an empty vector to store overlap measures for the current cell type
        bhattacharyya_coef <- numeric(length(sample_names))
        hellinger_dist <- numeric(length(sample_names))
        
        # Iterate over each sample
        for (i in 1:length(sample_names)) {
            
            # Extract distances from the current sample to reference samples
            sample_distances <- distance_data[[cell_type]]$query_to_ref_distances[sample_names[i], ]
            
            # Compute density of sample distances
            sample_density <- density(sample_distances)
            
            # Create a common grid for evaluating densities
            common_grid <- seq(min(min(ref_density$x), min(sample_density$x), 0), 
                               max(max(ref_density$x), max(sample_density$x)), length.out = 1000)
            
            # Interpolate densities onto the common grid
            ref_density_interp <- approxfun(ref_density$x, ref_density$y)(common_grid)
            ref_density_interp[is.na(ref_density_interp)] <- 0
            sample_density_interp <- approxfun(sample_density$x, sample_density$y)(common_grid)
            sample_density_interp[is.na(sample_density_interp)] <- 0
            
            # Compute and store Bhattacharyya coefficient/Hellinger distance
            bhattacharyya_coef[i] <- sum(sqrt(ref_density_interp * sample_density_interp) * mean(diff(common_grid)))
            hellinger_dist[i] <- sqrt(1 - sum(sqrt(ref_density_interp * sample_density_interp)) * mean(diff(common_grid)))
        }
        
        # Store overlap measures for the current cell type
        bhattacharyya_list[[cell_type]] <- bhattacharyya_coef
        hellinger_list[[cell_type]] <- hellinger_dist
    }
    
    # Return list with overlap measures
    bhattacharyya_coef <- data.frame(Sample = sample_names, bhattacharyya_list)
    hellinger_dist <- data.frame(Sample = sample_names, hellinger_list)
    return(list(bhattacharyya_coef = bhattacharyya_coef, 
                hellinger_dist = hellinger_dist))
}


