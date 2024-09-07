#' @title Calculate Nearest Neighbor Diagnostics for Cell Type Classification
#'
#' @description 
#' This function computes the probabilities for each query cell of belonging to either the reference or query dataset for 
#' each cell type using nearest neighbor analysis.
#'
#' @details 
#' The function conducts PCA on both the query and reference datasets to reduce dimensionality. It then compares each query 
#' cell to its nearest neighbors in the reference dataset to estimate the probability of membership in each cell type. Sample sizes 
#' between datasets are balanced using data augmentation if necessary.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A vector specifying the subset of principal components to use in the analysis. Default is 1:5
#' @param n_neighbor An integer specifying the number of nearest neighbors to consider. Default is 20.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#'
#' @return A list where each element corresponds to a cell type and contains:
#' \item{n_neighbor}{The number of nearest neighbors considered.}
#' \item{n_query}{The number of cells in the query dataset for each cell type.}
#' \item{query_prob}{The average probability of each query cell belonging to the reference dataset.}
#' The list is assigned the class \code{"calculateNearestNeighborProbabilities"}.
#' Each element in the list is named after the respective cell type.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.calculateNearestNeighborProbabilities}}
#' 
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#' 
#' # Project the query data onto PCA space of reference
#' nn_output <- calculateNearestNeighborProbabilities(query_data = query_data, 
#'                                                    reference_data = reference_data,
#'                                                    query_cell_type_col = "SingleR_annotation", 
#'                                                    ref_cell_type_col = "expert_annotation",
#'                                                    pc_subset = 1:10,
#'                                                    n_neighbor = 20)
#' 
#' # Plot output
#' plot(nn_output, cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"))
#' 
#' @importFrom stats dnorm
#' 
# Function to get probabilities for each cell of belonging to reference or query dataset for each cell type
calculateNearestNeighborProbabilities <- function(query_data, 
                                                  reference_data,
                                                  query_cell_type_col, 
                                                  ref_cell_type_col,
                                                  cell_types = NULL,
                                                  pc_subset = 1:5,
                                                  n_neighbor = 20,
                                                  assay_name = "logcounts"){
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)
    
    # Check if n_neighbor is a positive integer
    if (!is.numeric(n_neighbor) || n_neighbor <= 0 || 
        n_neighbor != as.integer(n_neighbor)) {
        stop("\'n_neighbor\' must be a positive integer.")
    }
    
    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }
    
    # Get PCA data
    pca_output <- projectPCA(query_data = query_data, 
                             reference_data = reference_data,
                             pc_subset = pc_subset,
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col,
                             assay_name = assay_name)
    
    # Initialize list to store probabilities
    probabilities <- vector("list", length = length(cell_types))
    names(probabilities) <- cell_types

    # Loop through each cell type
    for (cell_type in cell_types) {
        
        # Extract PCA-reduced data for the current cell type
        ref_pca_cell_type <- pca_output[which(
            pca_output[["dataset"]] == "Reference" & 
                pca_output[["cell_type"]] == cell_type), 
            paste0("PC", pc_subset)]
        query_pca_cell_type <- pca_output[which(
            pca_output[["dataset"]] == "Query" &
                pca_output[["cell_type"]] == cell_type), 
            paste0("PC", pc_subset)]
        
        # Combine reference and query data for the current cell type
        combined_data_cell_type <- rbind(ref_pca_cell_type, 
                                         query_pca_cell_type)
        
        # Number of cells for reference and query datasets
        n_ref <- nrow(ref_pca_cell_type)
        n_query <- nrow(query_pca_cell_type)
        
        # Data augmentation to balance sample size of datasets
        if(n_ref > n_query){
            
            combined_data_cell_type <- rbind(
                combined_data_cell_type,
                query_pca_cell_type[sample(seq_len(n_query), 
                                           n_ref - n_query, 
                                           replace = TRUE),])
            combined_data_cell_type <- data.frame(
                combined_data_cell_type, 
                data_type = rep(c("Reference", "Query"), each = n_ref))
        } else if (n_query > n_ref){
            
            combined_data_cell_type <- rbind(
                combined_data_cell_type,
                ref_pca_cell_type[sample(seq_len(n_ref), 
                                         n_query - n_ref, replace = TRUE),])
            combined_data_cell_type <- data.frame(
                combined_data_cell_type, 
                data_type = c(rep("Reference", n_ref), 
                              rep("Query", n_query),
                              rep("Reference", n_query - n_ref)))
        }
        
        # Perform nearest neighbors search
        dist_mat <- as.matrix(dist(
            combined_data_cell_type[, paste0("PC", pc_subset)]))
        neighbors_indices <- t(apply(dist_mat, 1, function(row, n_neighbor) {
            return(order(row)[2:(n_neighbor + 1)])},
        n_neighbor = n_neighbor))
        prob_query <- apply(neighbors_indices, 1, 
                            function(x, data_type) {
                                mean(data_type[x] == "Query")}, 
                            data_type = combined_data_cell_type$data_type)

        # Store the data for the cell type
        probabilities[[cell_type]] <- vector("list", length = 3)
        names(probabilities[[cell_type]]) <- c("n_neighbor", 
                                               "n_query", 
                                               "query_prob")
        probabilities[[cell_type]][["n_neighbor"]] <- n_neighbor
        probabilities[[cell_type]][["n_query"]] <- nrow(query_pca_cell_type)
        probabilities[[cell_type]][["query_prob"]] <- mean(prob_query)
    }
    
    # Creating class for output
    class(probabilities) <- c(class(probabilities), 
                              "calculateNearestNeighborProbabilities")
    
    # Return the list of probabilities
    return(probabilities)
}

