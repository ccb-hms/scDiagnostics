#' @title Calculate the Overlap Coefficient for Highly Variable Genes
#'
#' @description
#' Calculates the overlap coefficient between the sets of highly variable genes from a reference dataset and a query dataset.
#'
#' @details
#' The overlap coefficient measures the similarity between two gene sets, indicating how well-aligned
#' reference and query datasets are in terms of their highly variable genes. This metric is
#' useful in single-cell genomics to understand the correspondence between different datasets.
#'
#' The coefficient is calculated using the formula:
#'
#' \deqn{Coefficient(X, Y) = \frac{|X \cap Y|}{min(|X|, |Y|)}}
#'
#' where X and Y are the sets of highly variable genes from the reference and query datasets, respectively,
#' \eqn{|X \cap Y|} is the number of genes common to both \eqn{X} and \eqn{Y}, and \eqn{min(|X|, |Y|)} is the size of the
#' smaller set among \eqn{X} and \eqn{Y}.
#'
#' @param reference_genes A character vector of highly variable genes from the reference dataset.
#' @param query_genes A character vector of highly variable genes from the query dataset.
#'
#' @return Overlap coefficient, a value between 0 and 1, where 0 indicates no overlap
#' and 1 indicates complete overlap of highly variable genes between datasets.
#'
#' @export
#'
#' @references
#' Luecken et al. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19:41-50, 2022.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Selecting highly variable genes
#' ref_var <- scran::getTopHVGs(reference_data, n = 500)
#' query_var <- scran::getTopHVGs(query_data, n = 500)
#' overlap_coefficient <- calculateHVGOverlap(reference_genes = ref_var,
#'                                            query_genes = query_var)
#' overlap_coefficient
#'
# Function to calculate overlap between HVGs of reference and query datasets
calculateHVGOverlap <- function(reference_genes,
                                query_genes) {

    # Sanity checks
    if (!is.vector(reference_genes) || !is.character(reference_genes)) {
        stop("reference_genes must be a character vector.")
    }
    if (!is.vector(query_genes) || !is.character(query_genes)) {
        stop("query_genes must be a character vector.")
    }
    if (length(reference_genes) == 0 || length(query_genes) == 0) {
        stop("Input vectors must not be empty.")
    }

    # Calculate the intersection of highly variable genes
    common_genes <- intersect(reference_genes, query_genes)

    # Calculate the size of the intersection
    intersection_size <- length(common_genes)

    # Calculate the size of the smaller set
    min_size <- min(length(reference_genes), length(query_genes))

    # Compute the overlap coefficient
    overlap_coefficient <- intersection_size / min_size
    overlap_coefficient <- round(overlap_coefficient, 2)

    # Return the overlap coefficient
    return(overlap_coefficient)
}
