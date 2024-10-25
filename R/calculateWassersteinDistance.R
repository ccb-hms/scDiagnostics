#' @title Compute Wasserstein Distances Between Query and Reference Datasets
#'
#' @description
#' This function calculates Wasserstein distances between a query dataset and a reference dataset,
#' as well as within the reference dataset itself, after projecting them into a shared PCA space.
#'
#' @details
#' The function begins by projecting the query dataset onto the PCA space defined by the reference dataset.
#' It then computes Wasserstein distances between randomly sampled pairs within the reference dataset
#' to create a null distribution. Similarly, it calculates distances between the reference and query datasets.
#' The function assesses overall differences in distances to understand the variation between the datasets.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing a numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object with a numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies cell types.
#' @param pc_subset A numeric vector specifying which principal components to use. Default is \code{1:10}.
#' @param n_resamples An integer specifying the number of resamples to generate the null distribution. Default is \code{300}.
#' @param assay_name The name of the assay to use for computations. Default is \code{"logcounts"}.
#'
#' @return A list with the following components:
#' \item{null_dist}{A numeric vector of Wasserstein distances computed from resampled pairs within the reference dataset.}
#' \item{query_dist}{The mean Wasserstein distance between the query dataset and the reference dataset.}
#' \item{cell_type}{A character vector containing the unique cell types present in the reference dataset.}
#'
#' @references
#' Schuhmacher, D., Bernhard, S., & Book, M. (2019). "A Review of Approximate Transport in Machine Learning".
#' In \emph{Journal of Machine Learning Research} (Vol. 20, No. 117, pp. 1-61).
#'
#' @export
#'
#' @seealso \code{\link{plot.calculateWassersteinDistanceObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Extract CD4 cells
#' ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
#' query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]
#'
#' # Selecting highly variable genes (can be customized by the user)
#' ref_top_genes <- scran::getTopHVGs(ref_data_subset, n = 500)
#' query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_top_genes, query_top_genes)
#' ref_data_subset <- ref_data_subset[common_genes,]
#' query_data_subset <- query_data_subset[common_genes,]
#'
#' # Run PCA on reference data
#' ref_data_subset <- scater::runPCA(ref_data_subset)
#'
#' # Compute Wasserstein distances and compare using quantile-based permutation test
#' wasserstein_data <- calculateWassersteinDistance(query_data = query_data_subset,
#'                                                  reference_data = ref_data_subset,
#'                                                  query_cell_type_col = "expert_annotation",
#'                                                  ref_cell_type_col = "expert_annotation",
#'                                                  pc_subset = 1:5,
#'                                                  n_resamples = 100)
#' plot(wasserstein_data)
#'
#' @importFrom stats quantile
#'
# Function to generate density of Wasserstein distances under null distribution
calculateWassersteinDistance <- function(query_data,
                                         reference_data,
                                         ref_cell_type_col,
                                         query_cell_type_col,
                                         pc_subset = 1:5,
                                         n_resamples = 300,
                                         assay_name = "logcounts"){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  unique_cell_type = TRUE,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Check if n_resamples is a positive integer
    if (!inherits(n_resamples, "numeric")) {
        stop("\'n_resamples\' should be numeric.")
    } else if (any(!n_resamples == floor(n_resamples), n_resamples < 1)) {
        stop("\'n_resamples\' should be an integer, greater than zero.")
    }

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             pc_subset = pc_subset,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             assay_name = assay_name)

    # Get sample size for Wasserstein null distribution
    n_null <- min(floor(ncol(reference_data)/2), ncol(query_data), 500)

    # Extract variance explained
    weights <- attributes(reducedDim(
        reference_data, "PCA"))[["varExplained"]][pc_subset] /
        sum(attributes(reducedDim(
            reference_data, "PCA"))[["varExplained"]][pc_subset])

    # Compute reference-reference PCA weighted distances
    pca_ref <- pca_output[pca_output$dataset == "Reference",
                          paste0("PC", pc_subset)]
    pca_ref_weighted <- t(apply(pca_ref, 1,
                                function(x, weights) return(x * weights),
                                weights = sqrt(weights)))
    weighted_dist_ref <- as.matrix(dist(pca_ref_weighted))

    # Computing Wasserstein distances of null distribution
    null_dist <- numeric(n_resamples)
    prob_masses <- rep(1/n_null, n_null)
    for(iter in seq_len(n_resamples)){

        sample_ref_1 <- sample(seq_len(nrow(pca_ref)), n_null, replace = FALSE)
        sample_ref_2 <- sample(seq_len(nrow(pca_ref))[-sample_ref_1],
                               n_null, replace = FALSE)
        cost_mat <- weighted_dist_ref[sample_ref_1, sample_ref_2]
        opt_plan <- transport::transport(prob_masses, prob_masses,
                                         costm = cost_mat)
        null_dist[iter] <- transport::wasserstein(prob_masses,
                                                  prob_masses,
                                                  tplan = opt_plan,
                                                  costm = cost_mat)
    }

    # Compute reference-query PCA weighted distances
    pca_query <- pca_output[pca_output$dataset == "Query",
                            paste0("PC", pc_subset)]
    pca_query_weighted <- t(apply(pca_query, 1,
                                  function(x, weights) return(x * weights),
                                  weights = sqrt(weights)))
    weighted_dist_query <- outer(rowSums(pca_ref_weighted^2),
                                 rowSums(pca_query_weighted^2), "+") -
        2 * pca_ref_weighted %*% t(pca_query_weighted)

    # Computing Wasserstein distances for query data
    query_dist <- numeric(n_resamples)
    for(iter in seq_len(n_resamples)){

        sample_ref <- sample(seq_len(nrow(pca_ref)),
                             n_null, replace = FALSE)
        sample_query <- sample(seq_len(nrow(pca_query)),
                               n_null, replace = FALSE)
        cost_mat <- weighted_dist_query[sample_ref, sample_query]
        opt_plan <- transport::transport(prob_masses,
                                         prob_masses,
                                         costm = cost_mat)
        query_dist[iter] <- transport::wasserstein(prob_masses,
                                                   prob_masses,
                                                   tplan = opt_plan,
                                                   costm = cost_mat)
    }

    # Return the results
    wasserstein_data <- list(
        null_dist = null_dist,
        query_dist = mean(query_dist),
        cell_type = unique(reference_data[[ref_cell_type_col]])
    )
    class(wasserstein_data) <- c(class(wasserstein_data),
                                 "calculateWassersteinDistanceObject")
    return(wasserstein_data)
}







