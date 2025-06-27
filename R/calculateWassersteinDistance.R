#' @title Compute Wasserstein Distance Distributions Between Query and Reference Datasets
#'
#' @description
#' This function calculates distributions of Wasserstein distances between reference-reference
#' pairs and reference-query pairs for each specified cell type, after projecting them into
#' a shared PCA space. It then computes the probability of superiority to assess whether
#' reference-query distances tend to be larger than reference-reference distances.
#'
#' @details
#' The function projects the query dataset onto the PCA space defined by the reference dataset.
#' For each cell type, it computes two distributions: (1) Wasserstein distances between randomly
#' sampled pairs within the reference dataset (null distribution), and (2) Wasserstein distances
#' between reference and query dataset samples. It then calculates the probability of superiority,
#' which represents the probability that a randomly selected ref-query distance is larger than
#' a randomly selected ref-ref distance.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing a numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object with a numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies cell types.
#' @param cell_types A character vector specifying the cell types to include in the analysis. If NULL, all common cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to use. Default is \code{1:10}.
#' @param n_resamples An integer specifying the number of resamples to generate each distribution. Default is \code{300}.
#' @param assay_name The name of the assay to use for computations. Default is \code{"logcounts"}.
#'
#' @return A list with the following components:
#' \item{ref_ref_dist}{A named list of numeric vectors containing Wasserstein distances computed from resampled pairs within the reference dataset for each cell type.}
#' \item{ref_query_dist}{A named list of numeric vectors containing Wasserstein distances between reference and query datasets for each cell type.}
#' \item{probability_superiority}{A named numeric vector showing the probability that ref-query distances are larger than ref-ref distances for each cell type.}
#' \item{cell_types}{A character vector containing the cell types analyzed.}
#'
#' @references
#' Schuhmacher, D., Bernhard, S., & Book, M. (2019). "A Review of Approximate Transport in Machine Learning".
#' In \emph{Journal of Machine Learning Research} (Vol. 20, No. 117, pp. 1-61).
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateWassersteinDistanceObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compute Wasserstein distance distributions for all cell types
#' wasserstein_data <- calculateWassersteinDistance(query_data = query_data,
#'                                                  reference_data = reference_data,
#'                                                  query_cell_type_col = "expert_annotation",
#'                                                  ref_cell_type_col = "expert_annotation",
#'                                                  pc_subset = 1:10,
#'                                                  n_resamples = 100)
#' plot(wasserstein_data)
#'
#' @importFrom stats quantile
#'
calculateWassersteinDistance <- function(query_data,
                                         reference_data,
                                         ref_cell_type_col,
                                         query_cell_type_col,
                                         cell_types = NULL,
                                         pc_subset = 1:10,
                                         n_resamples = 300,
                                         assay_name = "logcounts"){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Check if n_resamples is a positive integer
    if (!inherits(n_resamples, "numeric")) {
        stop("\'n_resamples\' should be numeric.")
    } else if (any(!n_resamples == floor(n_resamples), n_resamples < 1)) {
        stop("\'n_resamples\' should be an integer, greater than zero.")
    }

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             pc_subset = pc_subset,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             assay_name = assay_name)

    # Split by cell type
    cell_list <- split(pca_output, pca_output[["cell_type"]])

    # Extract variance explained for weighting
    weights <- attributes(reducedDim(
        reference_data, "PCA"))[["varExplained"]][pc_subset] /
        sum(attributes(reducedDim(
            reference_data, "PCA"))[["varExplained"]][pc_subset])

    # Initialize results
    ref_ref_dist <- list()
    ref_query_dist <- list()
    probability_superiority <- numeric(length(cell_types))
    names(probability_superiority) <- cell_types

    # Process each cell type
    for (cell_type in cell_types) {

        # Skip if cell type not present in data
        if (!cell_type %in% names(cell_list)) {
            warning(paste("Cell type", cell_type, "not found in data. Skipping."))
            next
        }

        cell_data <- cell_list[[cell_type]]
        ref_indices <- cell_data[["dataset"]] == "Reference"
        query_indices <- cell_data[["dataset"]] == "Query"

        # Check if we have both reference and query cells for this cell type
        if (sum(ref_indices) == 0 || sum(query_indices) == 0) {
            warning(paste("Cell type", cell_type, "missing in reference or query data. Skipping."))
            next
        }

        # Get sample size for Wasserstein distributions
        n_sample <- min(floor(sum(ref_indices)/2), sum(query_indices), 200)

        if (n_sample < 10) {
            warning(paste("Too few cells for cell type", cell_type, ". Skipping."))
            next
        }

        # Extract PCA data for this cell type
        pca_ref <- as.matrix(cell_data[ref_indices, paste0("PC", pc_subset)])
        pca_query <- as.matrix(cell_data[query_indices, paste0("PC", pc_subset)])

        # Apply variance weighting
        pca_ref_weighted <- t(apply(pca_ref, 1,
                                    function(x, weights) return(x * weights),
                                    weights = sqrt(weights)))
        pca_query_weighted <- t(apply(pca_query, 1,
                                      function(x, weights) return(x * weights),
                                      weights = sqrt(weights)))

        # Compute reference-reference weighted distances (full distance matrix)
        weighted_dist_ref <- as.matrix(dist(pca_ref_weighted))

        # Compute reference-query weighted distances (full distance matrix)
        weighted_dist_query <- outer(rowSums(pca_ref_weighted^2),
                                     rowSums(pca_query_weighted^2), "+") -
            2 * pca_ref_weighted %*% t(pca_query_weighted)

        # Computing reference-reference Wasserstein distance distribution
        ref_ref_distances <- numeric(n_resamples)
        prob_masses <- rep(1/n_sample, n_sample)

        for(iter in seq_len(n_resamples)){
            sample_ref_1 <- sample(seq_len(nrow(pca_ref)), n_sample, replace = FALSE)
            sample_ref_2 <- sample(seq_len(nrow(pca_ref))[-sample_ref_1],
                                   n_sample, replace = FALSE)
            cost_mat <- weighted_dist_ref[sample_ref_1, sample_ref_2]
            opt_plan <- transport::transport(prob_masses, prob_masses,
                                             costm = cost_mat)
            ref_ref_distances[iter] <- transport::wasserstein(prob_masses,
                                                              prob_masses,
                                                              tplan = opt_plan,
                                                              costm = cost_mat)
        }

        # Computing reference-query Wasserstein distance distribution
        ref_query_distances <- numeric(n_resamples)

        for(iter in seq_len(n_resamples)){
            sample_ref <- sample(seq_len(nrow(pca_ref)), n_sample, replace = FALSE)
            sample_query <- sample(seq_len(nrow(pca_query)), n_sample, replace = FALSE)
            cost_mat <- weighted_dist_query[sample_ref, sample_query]
            opt_plan <- transport::transport(prob_masses, prob_masses,
                                             costm = cost_mat)
            ref_query_distances[iter] <- transport::wasserstein(prob_masses,
                                                                prob_masses,
                                                                tplan = opt_plan,
                                                                costm = cost_mat)
        }

        # Store distributions for this cell type
        ref_ref_dist[[cell_type]] <- ref_ref_distances
        ref_query_dist[[cell_type]] <- ref_query_distances

        # Calculate probability of superiority
        # P(ref_query > ref_ref) when sampling one value from each distribution
        n_comparisons <- 0
        n_superiority <- 0

        for (rq_val in ref_query_distances) {
            for (rr_val in ref_ref_distances) {
                n_comparisons <- n_comparisons + 1
                if (rq_val > rr_val) {
                    n_superiority <- n_superiority + 1
                }
            }
        }

        probability_superiority[cell_type] <- n_superiority / n_comparisons
    }

    # Filter out cell types that were skipped
    processed_cell_types <- names(ref_ref_dist)
    probability_superiority <- probability_superiority[processed_cell_types]

    # Return the results
    wasserstein_data <- list(
        ref_ref_dist = ref_ref_dist,
        ref_query_dist = ref_query_dist,
        probability_superiority = probability_superiority,
        cell_types = processed_cell_types
    )
    class(wasserstein_data) <- c(class(wasserstein_data),
                                 "calculateWassersteinDistanceObject")
    return(wasserstein_data)
}
