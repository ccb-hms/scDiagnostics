#' Calculate the Overlap Coefficient for Highly Variable Genes
#'
#' @description This function calculates the overlap coefficient between the sets of highly variable genes 
#' from a reference dataset and a query dataset. The overlap coefficient measures the 
#' degree of overlap or similarity between these two sets of genes, pointing towards how well reference
#' and query datasets are aligned.
#'
#' @param reference_genes A vector of highly variable genes from the reference dataset.
#' @param query_genes A vector of highly variable genes from the query dataset.
#'
#' 
#' @return Overlap coefficient, a value between 0 and 1, where 0 indicates no overlap 
#'         and 1 indicates complete overlap of highly variable genes between datasets.
#' 
#' @examples
#' 
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(SingleR)
#' 
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#' 
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#' 
#' # log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#' 
#' # Selcting highly variable genes
#' 
#' ref_var <- getTopHVGs(ref_data, n=2000)
#' query_var <- getTopHVGs(query_data, n=2000)
#' 
#' overlap_coefficient <- calculateHVGOverlap(reference_genes = ref_var, 
#'                                           query_genes = query_var)
#' @export                                       
calculateHVGOverlap <- function(reference_genes, query_genes) {
  
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