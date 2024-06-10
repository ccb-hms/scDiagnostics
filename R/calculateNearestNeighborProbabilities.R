#' @title Calculate Nearest Neighbor Diagnostics for Cell Type Classification
#'
#' @description 
#' This function computes the probabilities for each sample of belonging to either the reference or query dataset for 
#' each cell type using nearest neighbor analysis.

#'
#' @details 
#' This function performs a nearest neighbor search to calculate the probability of each sample in the query dataset 
#' belonging to the reference dataset for each cell type. It uses principal component analysis (PCA) to reduce the dimensionality 
#' of the data before performing the nearest neighbor search. The function balances the sample sizes between the reference and query 
#' datasets by data augmentation if necessary.

#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param n_neighbor An integer specifying the number of nearest neighbors to consider. Default is 15.
#' @param n_components An integer specifying the number of principal components to use for dimensionality reduction. Default is 10.
#' @param pc_subset A vector specifying the subset of principal components to use in the analysis. Default is c(1:10).
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#'
#' @return A list where each element corresponds to a cell type and contains two vectors:
#' \item{prob_ref}{The probabilities of each query sample belonging to the reference dataset.}
#' \item{prob_query}{The probabilities of each query sample belonging to the query dataset.}
#' The list is assigned the class \code{"nearestNeighbotDiagnostics"}.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.calculateNearestNeighborProbabilities}}
#' 
#' @examples
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
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
#' ref_var <- getTopHVGs(ref_data, n = 500)
#' query_var <- getTopHVGs(query_data, n = 500)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#' 
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' # Project the query data onto PCA space of reference
#' nn_output <- calculateNearestNeighborProbabilities(query_data_subset, ref_data_subset,
#'                                         n_neighbor = 15, 
#'                                         n_components = 10,
#'                                         pc_subset = c(1:10),
#'                                         query_cell_type_col = "labels", 
#'                                         ref_cell_type_col = "reclustered.broad")
#' 
#' # Plot output
#' plot(nn_output, cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'      prob_type = "query")
#' 
#' 
# Function to get probabilities for each sample of belonging to reference or query dataset for each cell type
calculateNearestNeighborProbabilities <- function(query_data, reference_data,
                                                  n_neighbor = 15,
                                                  n_components = 10,
                                                  pc_subset = c(1:10),
                                                  query_cell_type_col, 
                                                  ref_cell_type_col){
    
    # Check if n_components is a positive integer
    if (!inherits(n_components, "numeric")) {
        stop("n_components should be numeric")
    } else if (any(!n_components == floor(n_components), n_components < 1)) {
        stop("n_components should be an integer, greater than zero.")
    }
    
    # Get PCA data
    pca_output <- projectPCA(query_data = query_data, reference_data = reference_data,
                             n_components = n_components,
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col, 
                             return_value = c("data.frame", "list")[2])
    
    # Initialize list to store probabilities
    probabilities <- list()
    
    # Get unique cell types
    cell_types <- na.omit(intersect(unique(query_data[[query_cell_type_col]]), 
                                    unique(reference_data[[ref_cell_type_col]])))

    # Loop through each cell type
    for (cell_type in cell_types) {
        
        # Extract PCA-reduced data for the current cell type
        ref_pca_cell_type <- pca_output$ref[which(reference_data[[ref_cell_type_col]] == cell_type), paste0("PC", pc_subset)]
        query_pca_cell_type <- pca_output$query[which(query_data[[query_cell_type_col]] == cell_type), paste0("PC", pc_subset)]
        
        # Combine reference and query data for the current cell type
        combined_data_cell_type <- rbind(ref_pca_cell_type, query_pca_cell_type)
        
        # Number of samples for reference and query datasets
        n_ref <- nrow(ref_pca_cell_type)
        n_query <- nrow(query_pca_cell_type)
        
        # Data augmentation to balance sample size of datasets
        if(n_ref > n_query){
            
            combined_data_cell_type <- rbind(combined_data_cell_type,
                                             query_pca_cell_type[sample(1:n_query, n_ref - n_query, replace = TRUE),])
            combined_data_cell_type <- data.frame(combined_data_cell_type, data_type = rep(c("Reference", "Query"), each = n_ref))
        } else if (n_query > n_ref){
            
            combined_data_cell_type <- rbind(combined_data_cell_type,
                                             ref_pca_cell_type[sample(1:n_ref, n_query - n_ref, replace = TRUE),])
            combined_data_cell_type <- data.frame(combined_data_cell_type, data_type = c(rep("Reference", n_ref), rep("Query", n_query),
                                                                                         rep("Reference", n_query - n_ref)))
        }
        
        # Perform nearest neighbors search
        dist_mat <- as.matrix(dist(combined_data_cell_type[, paste0("PC", pc_subset)]))
        neighbors_indices <- t(apply(dist_mat, 1, get_lowest_indices <- function(row, n_neighbor) {
            return(order(row)[2:(n_neighbor + 1)])
        }, n_neighbor = n_neighbor))
        prob_ref <- apply(neighbors_indices, 1, function(x, data_type) {mean(data_type[x] == "Reference")}, 
                          data_type = combined_data_cell_type$data_type)
        
        # Store the probabilities
        probabilities[[cell_type]] <- list()
        probabilities[[cell_type]]$prob_ref <- prob_ref
        probabilities[[cell_type]]$prob_query <- 1 - prob_ref
    }
    
    # Creating class for output
    class(probabilities) <- c(class(probabilities), "calculateNearestNeighborProbabilities")
    
    # Return the list of probabilities
    return(probabilities)
}

