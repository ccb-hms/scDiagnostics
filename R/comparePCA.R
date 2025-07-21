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
#' pair of PCs using vectorized operations for improved performance. The resulting matrix contains the similarity values,
#' where rows represent reference PCs and columns represent query PCs.
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
#' @param n_permutations Number of permutations for statistical significance testing. If 0, no permutation test is performed. Default is 0.
#'
#' @return A list containing:
#'   \item{similarity_matrix}{A matrix comparing the principal components of the reference and query datasets.}
#'   \item{top_variables}{A list containing the top loading variables for each PC pair comparison.}
#'   \item{p_values}{A matrix of permutation p-values (if n_permutations > 0).}
#'   \item{metric}{The similarity metric used.}
#'   \item{n_top_vars}{Number of top variables used.}
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
#'                              correlation_method = c("spearman", "pearson")[1],
#'                              n_permutation = 100)
#'
#' # Create the heatmap
#' plot(similarity_mat, show_significance = TRUE)
#'
# Function to compare PCA of the SCE objects
comparePCA <- function(reference_data,
                       query_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       pc_subset = 1:5,
                       n_top_vars = 50,
                       metric = c("cosine", "correlation"),
                       correlation_method = c("spearman", "pearson"),
                       n_permutations = 0) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  unique_cell_type = TRUE,
                  pc_subset_query = pc_subset,
                  pc_subset_ref = pc_subset,
                  common_rotation_genes = TRUE)

    # Match arguments
    metric <- match.arg(metric)
    correlation_method <- match.arg(correlation_method)

    # Enhanced input validation
    if (!is.numeric(n_top_vars) || n_top_vars <= 0 ||
        n_top_vars != as.integer(n_top_vars)) {
        stop("'n_top_vars' must be a positive integer.")
    }

    if (!is.numeric(n_permutations) || n_permutations < 0 ||
        n_permutations != as.integer(n_permutations)) {
        stop("'n_permutations' must be a non-negative integer.")
    }

    # Extract PCA data from reference and query data
    ref_rotation <-
        attributes(reducedDim(reference_data, "PCA"))[["rotation"]][, pc_subset]
    query_rotation <-
        attributes(reducedDim(query_data, "PCA"))[["rotation"]][, pc_subset]
    query_rotation <-
        query_rotation[match(rownames(ref_rotation), rownames(query_rotation)), ]

    # Check if n_top_vars is reasonable given the number of genes
    n_genes <- nrow(ref_rotation)
    if (n_top_vars > n_genes) {
        warning(sprintf("n_top_vars (%d) is greater than the number of genes (%d). Using all genes.",
                        n_top_vars, n_genes))
        n_top_vars <- n_genes
    }

    # Function to identify high-loading variables for each PC (vectorized)
    .getHighLoadingVars <- function(rotation_mat, n_top_vars) {
        abs_loadings <- abs(rotation_mat)
        high_loading_vars <- apply(abs_loadings, 2, function(pc_loadings) {
            names(sort(pc_loadings, decreasing = TRUE))[seq_len(n_top_vars)]
        })
        return(high_loading_vars)
    }

    # Get top variables for each PC
    top_ref <- .getHighLoadingVars(ref_rotation, n_top_vars)
    top_query <- .getHighLoadingVars(query_rotation, n_top_vars)

    # Store top variables information for each comparison
    top_variables_info <- list()

    # Initialize similarity matrix
    similarity_matrix <- matrix(NA,
                                nrow = length(pc_subset),
                                ncol = length(pc_subset))

    # Vectorized similarity computation function
    .computeSimilarityMatrix <- function(ref_rot, query_rot, top_ref, top_query,
                                         metric, correlation_method) {
        sim_mat <- matrix(NA, nrow = ncol(ref_rot), ncol = ncol(query_rot))
        top_vars_list <- list()

        for (i in seq_len(ncol(ref_rot))) {
            for (j in seq_len(ncol(query_rot))) {
                # Get union of top variables for this PC pair
                combination_union <- union(top_ref[, i], top_query[, j])
                top_vars_list[[paste(i, j, sep = "_")]] <- combination_union

                # Extract relevant loadings
                ref_loadings <- ref_rot[combination_union, i]
                query_loadings <- query_rot[combination_union, j]

                # Compute similarity based on metric
                if (metric == "cosine") {
                    # Vectorized cosine similarity
                    dot_product <- sum(ref_loadings * query_loadings)
                    norm_ref <- sqrt(sum(ref_loadings^2))
                    norm_query <- sqrt(sum(query_loadings^2))
                    sim_mat[i, j] <- dot_product / (norm_ref * norm_query)
                } else if (metric == "correlation") {
                    sim_mat[i, j] <- cor(ref_loadings, query_loadings,
                                         method = correlation_method)
                }
            }
        }

        return(list(similarity_matrix = sim_mat, top_variables = top_vars_list))
    }

    # Compute similarity matrix
    result <- .computeSimilarityMatrix(ref_rotation, query_rotation,
                                       top_ref, top_query,
                                       metric, correlation_method)
    similarity_matrix <- result[["similarity_matrix"]]
    top_variables_info <- result[["top_variables"]]

    # Perform permutation test if requested
    p_values <- NULL
    if (n_permutations > 0) {
        message(sprintf("Performing %d permutations for significance testing...",
                        n_permutations))

        # Store original similarities for comparison
        original_similarities <- similarity_matrix
        permuted_similarities <- array(NA, dim = c(nrow(similarity_matrix),
                                                   ncol(similarity_matrix),
                                                   n_permutations))

        for (perm in seq_len(n_permutations)) {
            # Permute gene labels
            permuted_genes <- sample(rownames(ref_rotation))
            ref_rotation_perm <- ref_rotation
            rownames(ref_rotation_perm) <- permuted_genes

            # Recompute similarities with permuted data
            perm_result <- .computeSimilarityMatrix(ref_rotation_perm, query_rotation,
                                                    top_ref, top_query,
                                                    metric, correlation_method)
            permuted_similarities[, , perm] <- perm_result[["similarity_matrix"]]
        }

        # Calculate p-values
        p_values <- matrix(NA, nrow = nrow(similarity_matrix),
                           ncol = ncol(similarity_matrix))
        for (i in seq_len(nrow(similarity_matrix))) {
            for (j in seq_len(ncol(similarity_matrix))) {
                observed_sim <- abs(original_similarities[i, j])
                permuted_sims <- abs(permuted_similarities[i, j, ])
                p_values[i, j] <- (sum(permuted_sims >= observed_sim) + 1) /
                    (n_permutations + 1)
            }
        }
    }

    # Add informative row and column names
    ref_var_explained <-
        attributes(reducedDim(reference_data, "PCA"))[["percentVar"]][pc_subset]
    query_var_explained <-
        attributes(reducedDim(query_data, "PCA"))[["percentVar"]][pc_subset]

    rownames(similarity_matrix) <-
        paste0("Ref PC", pc_subset, " (", round(ref_var_explained, 1), "%)")
    colnames(similarity_matrix) <-
        paste0("Query PC", pc_subset, " (", round(query_var_explained, 1), "%)")

    if (!is.null(p_values)) {
        rownames(p_values) <- rownames(similarity_matrix)
        colnames(p_values) <- colnames(similarity_matrix)
    }

    # Create comprehensive output object
    output <- list(
        similarity_matrix = similarity_matrix,
        top_variables = top_variables_info,
        p_values = p_values,
        metric = metric,
        correlation_method = if(metric == "correlation") correlation_method else NULL,
        n_top_vars = n_top_vars,
        n_permutations = n_permutations,
        pc_subset = pc_subset
    )

    # Update class
    class(output) <- c("comparePCAObject", class(output))

    return(output)
}
