#' @title Compare Subspaces Spanned by Top Principal Components
#'
#' @description
#' This function compares the subspace spanned by the top principal components (PCs) in a reference dataset to that
#' in a query dataset. It computes the cosine similarity between the loadings of the top variables for each PC in
#' both datasets and provides a weighted cosine similarity score.
#'
#' @details
#' This function compares the subspace spanned by the top principal components (PCs) in a reference dataset
#' to that in a query dataset. It first computes the cosine similarity between the loadings of the top variables
#' for each PC in both datasets. The top cosine similarity scores are then selected, and their corresponding PC
#' indices are stored. Additionally, the function calculates the average percentage of variance explained by the
#' selected top PCs. Finally, it computes a weighted cosine similarity score based on the top cosine similarities
#' and the average percentage of variance explained.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) to compare. Default is the first five PCs.
#' @param n_top_vars An integer indicating the number of top loading variables to consider for each PC. Default is 50.
#'
#' @return A list containing the following components:
#'   \item{cosine_similarity}{A numeric vector of cosine values of principal angles.}
#'   \item{cosine_id}{A matrix showing which reference and query PCs were matched.}
#'   \item{var_explained_ref}{A numeric vector of variance explained by reference PCs.}
#'   \item{var_explained_query}{A numeric vector of variance explained by query PCs.}
#'   \item{var_explained_avg}{A numeric vector of average variance explained by each PC pair.}
#'   \item{weighted_cosine_similarity}{A numeric value representing the weighted cosine similarity.}
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.comparePCASubspaceObject}}
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
#' # Compare PCA subspaces
#' subspace_comparison <- comparePCASubspace(query_data = query_data_subset,
#'                                           reference_data = ref_data_subset,
#'                                           query_cell_type_col = "expert_annotation",
#'                                           ref_cell_type_col = "expert_annotation",
#'                                           n_top_vars = 50,
#'                                           pc_subset = 1:5)
#'
#' # Plot output for PCA subspace comparison
#' plot(subspace_comparison)
#'
# Function to compare subspace spanned by top PCs in reference and query datasets
comparePCASubspace <- function(reference_data,
                               query_data,
                               query_cell_type_col,
                               ref_cell_type_col,
                               pc_subset = 1:5,
                               n_top_vars = 50){

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

    # Compute the cosine similarity (cosine of principal angle)
    cosine_similarity_result <- comparePCA(query_data = query_data,
                                           reference_data = reference_data,
                                           query_cell_type_col = query_cell_type_col,
                                           ref_cell_type_col = ref_cell_type_col,
                                           pc_subset = pc_subset,
                                           n_top_vars = n_top_vars,
                                           metric = "cosine")

    # Extract similarity matrix from the result
    cosine_similarity <- cosine_similarity_result[["similarity_matrix"]]

    # Vector to store top cosine similarities
    top_cosine <- numeric(length(pc_subset))
    # Matrix to store PC IDs for each top cosine similarity
    cosine_id <- matrix(NA, nrow = length(pc_subset), ncol = 2)
    colnames(cosine_id) <- c("Ref", "Query")

    # Looping to store top cosine similarities and PC IDs
    for(id in seq_len(length(pc_subset))){

        # Store data for top cosine
        top_ref <- which.max(apply(abs(cosine_similarity), 1, max))
        top_query <- which.max(abs(cosine_similarity)[top_ref,])
        top_cosine[id] <- abs(cosine_similarity)[top_ref, top_query]
        cosine_id[id,] <- c(top_ref, top_query)

        # Remove as candidate
        cosine_similarity[top_ref,] <- 0
        cosine_similarity[, top_query] <- 0
    }

    # Vector of variance explained - FIXED BUG AND ADDED INDIVIDUAL VALUES
    var_explained_ref_all <- attributes(
        reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset]
    var_explained_query_all <- attributes(
        reducedDim(query_data, "PCA"))[["percentVar"]][pc_subset]

    # Get variance explained for the matched PCs
    var_explained_ref <- var_explained_ref_all[cosine_id[, 1]]
    var_explained_query <- var_explained_query_all[cosine_id[, 2]]
    var_explained_avg <- (var_explained_ref + var_explained_query) / 2

    # Weighted cosine similarity score
    weighted_cosine_similarity <- sum(top_cosine * var_explained_avg)/100

    # Update class of return output
    output <- list(cosine_similarity = top_cosine,
                   cosine_id = cosine_id,
                   var_explained_ref = var_explained_ref,
                   var_explained_query = var_explained_query,
                   var_explained_avg = var_explained_avg,
                   weighted_cosine_similarity = weighted_cosine_similarity)
    class(output) <- c(class(output), "comparePCASubspaceObject")

    # Return cosine similarity output
    return(output)
}
