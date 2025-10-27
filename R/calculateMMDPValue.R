#' @title Calculate Maximum Mean Discrepancy P-Values for Two-Sample Comparison
#'
#' @description
#' This function performs the Maximum Mean Discrepancy (MMD) test for comparing
#' distributions between two samples in PCA space using a custom implementation
#' with permutation testing for better sensitivity.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Projects the data into the PCA space.
#'   \item Subsets the data to the specified cell types and principal components.
#'   \item Performs a custom MMD test with permutation-based p-values for each cell type.
#' }
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#' @param n_permutation Number of permutations for p-value calculation. Default is 100.
#' @param kernel_type Type of kernel to use. Options are "gaussian" (default) or "linear".
#' @param sigma Bandwidth parameter for Gaussian kernel. If NULL, uses median heuristic.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#'
#' @return A named vector of p-values from the MMD test for each cell type.
#'
#' @references Gretton, A., Borgwardt, K. M., Rasch, M. J., Schölkopf, B., & Smola, A. (2012).
#' "A kernel two-sample test". Journal of Machine Learning Research, 13(1), 723-773.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Calculate MMD p-values (with query data)
#' mmd_test <- calculateMMDPValue(reference_data = reference_data,
#'                               query_data = query_data,
#'                               ref_cell_type_col = "expert_annotation",
#'                               query_cell_type_col = "SingleR_annotation",
#'                               cell_types = c("CD4", "CD8"),
#'                               pc_subset = 1:5,
#'                               n_permutation = 30)
#' mmd_test
#'
#' @importFrom stats median
#'
# Function to perform MMD test for two-sample comparison
calculateMMDPValue <- function(query_data,
                               reference_data,
                               query_cell_type_col,
                               ref_cell_type_col,
                               cell_types = NULL,
                               pc_subset = 1:5,
                               n_permutation = 100,
                               kernel_type = "gaussian",
                               sigma = NULL,
                               assay_name = "logcounts",
                               max_cells_query = 5000,
                               max_cells_ref = 5000) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cell_type_col)
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = TRUE,
                                  n_cell_types = NULL)

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             cell_types = cell_types,
                             pc_subset = pc_subset,
                             assay_name = assay_name,
                             max_cells_ref = max_cells_ref,
                             max_cells_query = max_cells_query)

    # Set data for MMD test
    cell_list <- split(pca_output, pca_output[["cell_type"]])

    # Set the PC variables
    pc_vars <- paste0("PC", pc_subset)

    # Calculate MMD p-values for each cell type
    p_values <- vector("numeric", length = length(cell_types))
    names(p_values) <- cell_types

    for(cell_type in cell_types){
        dataset_ind <- cell_list[[cell_type]][, "dataset"] == "Reference"

        X <- as.matrix(cell_list[[cell_type]][dataset_ind, pc_vars])
        Y <- as.matrix(cell_list[[cell_type]][!dataset_ind, pc_vars])

        # Skip if either group is too small
        if(nrow(X) < 3 || nrow(Y) < 3) {
            p_values[cell_type] <- NA
            next
        }

        # Adaptive permutation testing - stop early if clear result
        observed_mmd <- computeMMDStatistic(X, Y, kernel_type, sigma)

        combined_data <- rbind(X, Y)
        n_total <- nrow(combined_data)
        n_X <- nrow(X)

        # Adaptive permutation with early stopping
        min_perms <- min(100, n_permutation)
        extreme_count <- 0

        for(i in seq_len(n_permutation)) {
            perm_indices <- sample.int(n_total)
            X_perm <-
                combined_data[perm_indices[seq_len(n_X)], , drop = FALSE]
            Y_perm <-
                combined_data[perm_indices[seq(n_X + 1, n_total)], , drop = FALSE]
            perm_stat <-
                computeMMDStatistic(X_perm, Y_perm, kernel_type, sigma)

            if(perm_stat >= observed_mmd) {
                extreme_count <- extreme_count + 1
            }

            # Early stopping if we have enough evidence
            if(i >= min_perms) {
                current_pval <- (extreme_count + 1) / (i + 1)
                # Stop early if p-value is clearly < 0.01 or > 0.1
                if(current_pval < 0.005 || current_pval > 0.15) {
                    p_values[cell_type] <- (extreme_count + 1) / (i + 1)
                    break
                }
            }

            # If we've done all permutations
            if(i == n_permutation) {
                p_values[cell_type] <- (extreme_count + 1) /
                    (n_permutation + 1)
            }
        }
    }

    return(p_values)
}

#' @title Compute Maximum Mean Discrepancy Statistic
#'
#' @description
#' Compute the Maximum Mean Discrepancy (MMD) statistic between two datasets
#' using either Gaussian or linear kernels for distribution comparison.
#'
#' @details
#' This function calculates the MMD statistic, which measures the distance between
#' two probability distributions by comparing their embeddings in a reproducing
#' kernel Hilbert space (RKHS). For the Gaussian kernel, an optimized median
#' heuristic is used to estimate the bandwidth parameter sigma when not provided.
#' The linear kernel provides a computationally faster alternative.
#'
#' The MMD statistic is computed as:
#' MMD^2 = E[k(X,X')] + E[k(Y,Y')] - 2*E[k(X,Y)]
#' where k is the chosen kernel function.
#'
#' @param X A numeric matrix representing the first dataset, where rows are observations
#' and columns are features.
#' @param Y A numeric matrix representing the second dataset, where rows are observations
#' and columns are features.
#' @param kernel_type A character string specifying the kernel type. Options are
#' "gaussian" for RBF kernel or "linear" for linear kernel. Default is "gaussian".
#' @param sigma A numeric value specifying the bandwidth parameter for the Gaussian kernel.
#' If NULL, it is estimated using the median heuristic. Default is NULL.
#'
#' @keywords internal
#'
#' @return A numeric value representing the MMD^2 statistic between the two datasets.
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to compute MMD statistic
computeMMDStatistic <- function(X, Y,
                                kernel_type = "gaussian",
                                sigma = NULL) {

    # Get sample sizes
    n <- nrow(X)
    m <- nrow(Y)

    # Gaussian kernel computation
    if (kernel_type == "gaussian") {

        # Estimate sigma using median heuristic if not provided
        if (is.null(sigma)) {
            # Efficient median heuristic using sampled subset of distances
            combined <- rbind(X, Y)
            n_combined <- nrow(combined)

            # Sample maximum of 1000 pairs to avoid computational burden
            max_pairs <- min(1000, n_combined * (n_combined - 1) / 2)
            sampled_dists <- numeric(max_pairs)

            for(k in seq_len(max_pairs)) {
                i <- sample.int(n_combined, 1)
                j <- sample.int(n_combined, 1)
                if(i != j) {
                    sampled_dists[k] <- sqrt(sum((combined[i, ] -
                                                      combined[j, ])^2))
                }
            }
            sigma <- median(sampled_dists[sampled_dists > 0])
        }

        # Pre-compute kernel parameter for efficiency
        sigma_sq_2_inv <- 1 / (2 * sigma^2)

        # Initialize kernel sums
        K_XX_sum <- 0
        K_YY_sum <- 0
        K_XY_sum <- 0

        # Compute K(X,X) sum - kernel evaluations within first dataset
        for(i in seq_len(n - 1)) {
            X_i <- X[i, ]
            for(j in seq(i + 1, n)) {
                dist_sq <- sum((X_i - X[j, ])^2)
                K_XX_sum <- K_XX_sum + 2 * exp(-dist_sq * sigma_sq_2_inv)
            }
        }

        # Compute K(Y,Y) sum - kernel evaluations within second dataset
        for(i in seq_len(m - 1)) {
            Y_i <- Y[i, ]
            for(j in seq(i + 1, m)) {
                dist_sq <- sum((Y_i - Y[j, ])^2)
                K_YY_sum <- K_YY_sum + 2 * exp(-dist_sq * sigma_sq_2_inv)
            }
        }

        # Compute K(X,Y) sum - kernel evaluations between datasets
        for(i in seq_len(n)) {
            X_i <- X[i, ]
            for(j in seq_len(m)) {
                dist_sq <- sum((X_i - Y[j, ])^2)
                K_XY_sum <- K_XY_sum + exp(-dist_sq * sigma_sq_2_inv)
            }
        }

    } else {
        # Linear kernel - computationally efficient alternative
        XX <- tcrossprod(X)
        YY <- tcrossprod(Y)
        XY <- tcrossprod(X, Y)

        # Compute kernel sums excluding diagonal elements for within-dataset terms
        K_XX_sum <- sum(XX) - sum(diag(XX))
        K_YY_sum <- sum(YY) - sum(diag(YY))
        K_XY_sum <- sum(XY)
    }

    # Calculate MMD^2 statistic using unbiased estimator
    mmd_stat <- K_XX_sum/(n*(n-1)) + K_YY_sum/(m*(m-1)) - 2*K_XY_sum/(n*m)

    return(mmd_stat)
}
