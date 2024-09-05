#' @title Compare Subspaces Spanned by Top Principal Components Using Canonical Correlation Analysis
#' 
#' @description 
#' This function compares the subspaces spanned by the top principal components (PCs) of the reference 
#' and query datasets using canonical correlation analysis (CCA). It calculates the canonical variables, 
#' correlations, and a similarity measure for the subspaces.
#'
#' @details
#' This function performs canonical correlation analysis (CCA) to compare the subspaces spanned by the 
#' top principal components (PCs) of the reference and query datasets. The function extracts the rotation 
#' matrices corresponding to the specified PCs and performs CCA on these matrices. It computes the canonical 
#' variables and their corresponding correlations. Additionally, it calculates a similarity measure for the 
#' canonical variables using cosine similarity. The output is a list containing the canonical coefficients 
#' for both datasets, the cosine similarity values, and the canonical correlations.

#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) 
#' to compare. Default is the first five PCs.
#' @param n_top_vars An integer indicating the number of top loading variables to consider for each PC. Default is 25.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{coef_ref}{Canonical coefficients for the reference dataset.}
#'   \item{coef_query}{Canonical coefficients for the query dataset.}
#'   \item{cosine_similarity}{Cosine similarity values for the canonical variables.}
#'   \item{correlations}{Canonical correlations between the reference and query datasets.}
#' }
#'
#' @export
#' 
#' @references
#' Hotelling, H. (1936). "Relations between two sets of variates". *Biometrika*, 28(3/4), 321-377. doi:10.2307/2333955.
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.compareCCA}}
#' 
#' @examples
#' # Load libraries
#' library(scran)
#' library(scater)
#' 
#' # Load data
#' data("reference_data")
#' data("query_data")
#' 
#' # Extract CD4 cells
#' ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
#' query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_top_genes <- getTopHVGs(ref_data_subset, n = 500)
#' query_top_genes <- getTopHVGs(query_data_subset, n = 500)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_top_genes, query_top_genes)
#' ref_data_subset <- ref_data_subset[common_genes,]
#' query_data_subset <- query_data_subset[common_genes,]
#' 
#' # Run PCA on datasets separately
#' ref_data_subset <- runPCA(ref_data_subset)
#' query_data_subset <- runPCA(query_data_subset)
#' 
#' # Compare CCA
#' cca_comparison <- compareCCA(query_data = query_data_subset, 
#'                              reference_data = ref_data_subset, 
#'                              query_cell_type_col = "expert_annotation", 
#'                              ref_cell_type_col = "expert_annotation", 
#'                              pc_subset = 1:5, 
#'                              n_top_vars = 25)
#' 
#' # Visualize output of CCA comparison
#' plot(cca_comparison)
#' 
#' 
# Function to compare subspace spanned by top PCs in reference and query datasets
compareCCA <- function(query_data,
                       reference_data, 
                       query_cell_type_col, 
                       ref_cell_type_col, 
                       pc_subset = 1:5,
                       n_top_vars = 25){
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  unique_cell_type = TRUE,
                  pc_subset_query = pc_subset,
                  pc_subset_ref = pc_subset, 
                  common_rotation_genes = TRUE)
    
    # Check if n_top_vars is a positive integer
    if (!is.numeric(n_top_vars) || n_top_vars <= 0 || 
        n_top_vars != as.integer(n_top_vars)) {
        stop("\'n_top_vars\' must be a positive integer.")
    }
    
    # Extract the rotation matrices
    ref_rotation <- attributes(
        reducedDim(reference_data, "PCA"))[["rotation"]][, pc_subset]
    query_rotation <- attributes(
        reducedDim(query_data, "PCA"))[["rotation"]][, pc_subset]
    query_rotation <- query_rotation[match(rownames(ref_rotation), 
                                           rownames(query_rotation)), ]

    # Function to identify high-loading variables for each PC
    .getHighLoadingVars <- function(rotation_mat, n_top_vars) {
        high_loading_vars <- lapply(seq_len(ncol(rotation_mat)), function(pc) {
            abs_loadings <- abs(rotation_mat[, pc])
            top_vars <- names(sort(abs_loadings, 
                                   decreasing = TRUE))[seq_len(n_top_vars)]
            return(top_vars)
        })
        return(high_loading_vars)
    }
    
    # Get union of variables with highest loadings
    top_ref <- .getHighLoadingVars(ref_rotation, n_top_vars)
    top_query <- .getHighLoadingVars(query_rotation, n_top_vars)
    top_union <- unlist(lapply(seq_len(length(pc_subset)), 
                               function(i) return(union(top_ref[[i]], 
                                                        top_query[[i]]))))
    
    # Perform CCA
    cca_result <- cancor(ref_rotation, query_rotation)
    
    # Extract canonical variables and correlations
    canonical_ref <- cca_result$xcoef
    canonical_query <- cca_result$ycoef
    correlations <- cca_result$cor
    
    # Function to compute similarity measure (e.g., cosine similarity)
    .cosine_similarity <- function(u, v) {
        return(abs(sum(u * v)) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
    }
    
    # Compute similarities and account for correlations
    similarities <- rep(0, length(pc_subset))
    for (i in seq_len(length(pc_subset))) {
        similarities[i] <- .cosine_similarity(canonical_ref[, i], 
                                              canonical_query[, i])
    }
    
    # Update class of return output
    output <- list(coef_ref = canonical_ref,
                   coef_query = canonical_query,
                   cosine_similarity = similarities,
                   correlations = correlations)
    class(output) <- c(class(output), "compareCCA")

    # Return cosine similarity output
    return(output)
}

