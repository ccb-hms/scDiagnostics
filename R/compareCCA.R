#' @title Compare Canonical Correlation Analysis (CCA) between Query and Reference Data
#'
#' @description
#' This function performs Canonical Correlation Analysis (CCA) between two datasets (query and reference) after
#' performing PCA on each dataset. It projects the query data onto the PCA space of the reference data and then
#' computes the cosine similarity of the canonical correlation vectors between the two datasets.
#'
#' @details
#' The function performs the following steps:
#' 1. Projects the query data onto the PCA space of the reference data using the specified number of principal components.
#' 2. Downsamples the datasets to ensure an equal number of rows for CCA.
#' 3. Performs CCA on the projected datasets.
#' 4. Computes the cosine similarity between the canonical correlation vectors and extracts the canonical correlations.
#'
#' The cosine similarity provides a measure of alignment between the canonical correlation vectors of the two datasets.
#' Higher values indicate greater similarity.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs)
#' to compare. Default is the first five PCs.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
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
#' @seealso \code{\link{plot.compareCCAObject}}
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
#'                              pc_subset = 1:5)
#'
#' # Visualize output of CCA comparison
#' plot(cca_comparison)
#'
# Function to compare subspace spanned by top PCs in reference and query datasets
compareCCA <- function(query_data,
                       reference_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       pc_subset = 1:5,
                       assay_name = "logcounts"){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  unique_cell_type = TRUE,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             pc_subset = pc_subset,
                             assay_name = assay_name)

    # Get respective datasets
    ref_data <- pca_output[pca_output[["dataset"]] == "Reference", paste0("PC", pc_subset)]
    query_data <- pca_output[pca_output[["dataset"]] == "Query", paste0("PC", pc_subset)]

    # Downsample to ensure equal size of datasets
    n_samples <- min(nrow(ref_data), nrow(query_data))
    if (nrow(query_data) > n_samples) {
        query_data <- query_data[sample(nrow(query_data), n_samples), ]
    } else if (nrow(ref_data) > n_samples) {
        ref_data <- ref_data[sample(nrow(ref_data), n_samples), ]
    }

    # Perform CCA
    cca_result <- cancor(ref_data, query_data)

    # Extract canonical variables and correlations
    canonical_ref <- cca_result[["xcoef"]]
    canonical_query <- cca_result[["ycoef"]]
    correlations <- cca_result[["cor"]]

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
    class(output) <- c(class(output), "compareCCAObject")

    # Return cosine similarity output
    return(output)
}

