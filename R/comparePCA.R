#' @title Compare Principal Components Analysis (PCA) Results
#'
#' @description This function compares the principal components (PCs) obtained from separate PCA on reference and query
#' datasets for a single cell type using either cosine similarity or correlation.
#'
#' @details
#' This function compares the PCA results between the reference and query datasets by computing cosine
#' similarities or correlations between the loadings of top variables for each pair of principal components. It first
#' extracts the PCA rotation matrices from both datasets and identifies the top variables with highest loadings for
#' each PC. Then, it computes the cosine similarities or correlations between the loadings of top variables for each
#' pair of PCs. The resulting matrix contains the similarity values, where rows represent reference PCs and columns
#' represent query PCs.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) to compare. Default is the first five PCs.
#' @param n_top_vars An integer indicating the number of top loading variables to consider for each PC. Default is 50.
#' @param metric The similarity metric to use. It can be either "cosine" or "correlation". Default is "cosine".
#' @param correlation_method The correlation method to use if metric is "correlation". It can be "spearman"
#' or "pearson". Default is "spearman".
#'
#' @return A similarity matrix comparing the principal components of the reference and query datasets.
#' Each element (i, j) in the matrix represents the similarity between the i-th principal component
#' of the reference dataset and the j-th principal component of the query dataset.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.comparePCAObject}}
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
#' # Call the PCA comparison function
#' similarity_mat <- comparePCA(query_data = query_data_subset,
#'                              reference_data = ref_data_subset,
#'                              query_cell_type_col = "expert_annotation",
#'                              ref_cell_type_col = "expert_annotation",
#'                              pc_subset = 1:5,
#'                              n_top_vars = 50,
#'                              metric = c("cosine", "correlation")[1],
#'                              correlation_method = c("spearman", "pearson")[1])
#'
#' # Create the heatmap
#' plot(similarity_mat)
#'
# Compare PCA vectors of reference and query datasets for specific cell type.
comparePCA <- function(reference_data,
                       query_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       pc_subset = 1:5,
                       n_top_vars = 50,
                       metric = c("cosine", "correlation"),
                       correlation_method = c("spearman", "pearson")){

    # Match metric argument
    metric <- match.arg(metric)

    # Match correlation method argument
    correlation_method <- match.arg(correlation_method)

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

    # Extract PCA data from reference and query data
    ref_rotation <- attributes(
        reducedDim(reference_data, "PCA"))[["rotation"]][, pc_subset]
    query_rotation <- attributes(
        reducedDim(query_data, "PCA"))[["rotation"]][, pc_subset]
    query_rotation <- query_rotation[match(rownames(ref_rotation),
                                           rownames(query_rotation)),]

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
    top_union <- lapply(seq_len(length(pc_subset)),
                        function(i) return(union(top_ref[[i]], top_query[[i]])))

    # Initialize a matrix to store cosine similarities
    similarity_matrix <- matrix(NA, nrow = length(pc_subset),
                                ncol = length(pc_subset))

    if(metric == "cosine"){
        # Function to compute cosine similarity
        .cosine_similarity <- function(x, y) {
            sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
        }

        # Loop over each pair of columns and compute cosine similarity
        for (i in seq_len(length(pc_subset))) {
            for (j in seq_len(length(pc_subset))) {
                combination_union <- union(top_union[[i]], top_union[[j]])
                similarity_matrix[i, j] <- .cosine_similarity(
                    ref_rotation[combination_union, i],
                    query_rotation[combination_union, j])
            }
        }
    } else if(metric == "correlation"){
        # Loop over each pair of columns and compute cosine similarity
        for (i in seq_len(length(pc_subset))) {
            for (j in seq_len(length(pc_subset))) {
                combination_union <- union(top_union[[i]], top_union[[j]])
                similarity_matrix[i, j] <- cor(
                    ref_rotation[combination_union, i],
                    query_rotation[combination_union, j],
                    method = correlation_method)
            }
        }
    }

    # Add rownames and colnames with % of variance explained for each PC of each dataset
    rownames(similarity_matrix) <- paste0(
        "Ref PC", pc_subset, " (",
        round(attributes(reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset], 1), "%)")
    colnames(similarity_matrix) <- paste0(
        "Query PC", pc_subset, " (",
        round(attributes(reducedDim(query_data, "PCA"))[["percentVar"]][pc_subset], 1), "%)")

    # Update class of return output
    class(similarity_matrix) <- c(class(similarity_matrix), "comparePCAObject")

    # Return similarity matrix
    return(similarity_matrix)
}


