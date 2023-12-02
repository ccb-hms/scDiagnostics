#' Compute Average Pairwise Correlation between Cell Types
#'
#' Computes the average pairwise correlations between specified cell types 
#' in single-cell gene expression data.
#' 
#' @details This function operates on \code{\linkS4class{SingleCellExperiment}} objects, 
#' ideal for single-cell analysis workflows. It calculates pairwise correlations between query and 
#' reference cells using a specified correlation method, then averages these correlations for each 
#' cell type pair. This function aids in assessing the similarity between cells in reference and query datasets, 
#' providing insights into the reliability of cell type annotations in single-cell gene expression data.
#'
#' @param query_data  A \code{\linkS4class{SingleCellExperiment}} containing the single-cell 
#' expression data and metadata.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing the single-cell 
#' expression data and metadata.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to be analysed consider.
#' @param correlation_method The correlation method to use for calculating pairwise correlations.
#'
#' @return A matrix containing the average pairwise correlation values. 
#'         Rows and columns are labeled with the cell types. Each element 
#'         in the matrix represents the average correlation between a pair 
#'         of cell types.
#' 
#' @examples
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
#' # Get cell type scores using SingleR
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#'
#' # Compute Pairwise Correlations
#' # Note: The selection of highly variable genes and desired cell types may vary 
#' # based on user preference. 
#' # The cell type annotation method used in this example is SingleR. 
#' # User can use any other method for cell type annotation and provide 
#' # the corresponding labels in the metadata.
#'
#' # Selecting highly variable genes
#' ref_var <- getTopHVGs(ref_data, n = 2000)
#' query_var <- getTopHVGs(query_data, n = 2000)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#'
#' # Select desired cell types
#' selected_cell_types <- c("CD4", "CD8", "B_and_plasma")
#' ref_data_subset <- ref_data[common_genes, ref_data$reclustered.broad %in% selected_cell_types]
#' query_data_subset <- query_data[common_genes, query_data$reclustered.broad %in% selected_cell_types]
#'
#' # Compute pairwise correlations
#' cor_matrix_avg <- computeAveragePairwiseCorrelation(query_data = query_data_subset, 
#'                                                     reference_data = ref_data_subset, 
#'                                                     query_cell_type_col = "labels", 
#'                                                     ref_cell_type_col = "reclustered.broad", 
#'                                                     cell_types = selected_cell_types, 
#'                                                     correlation_method = "spearman")
#'
#' # Visualize the results using any visualization method of choice
#'
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor
#' @export
computeAveragePairwiseCorrelation <- function(query_data, 
                                              reference_data, 
                                              query_cell_type_col, 
                                              ref_cell_type_col, 
                                              cell_types, correlation_method) {
  # Sanity checks
  
  # Check if query_data is a SingleCellExperiment object
  if (!is(query_data, "SingleCellExperiment")) {
    stop("query_data must be a SingleCellExperiment object.")
  }
  
  # Check if reference_data is a SingleCellExperiment object
  if (!is(reference_data, "SingleCellExperiment")) {
    stop("reference_data must be a SingleCellExperiment object.")
  }
  
  # Check if query_cell_type_col is a valid column name in query_data
  if (!query_cell_type_col %in% colnames(colData(query_data))) {
    stop("query_cell_type_col: '", query_cell_type_col, "' is not a valid column name in query_data.")
  }
  
  # Check if ref_cell_type_col is a valid column name in reference_data
  if (!ref_cell_type_col %in% colnames(colData(reference_data))) {
    stop("ref_cell_type_col: '", ref_cell_type_col, "' is not a valid column name in reference_data.")
  }
  
  # Check if all cell_types are present in query_data
  if (!all(cell_types %in% unique(query_data[[query_cell_type_col]]))) {
    stop("One or more cell_types specified are not present in query_data.")
  }
  
  # Check if all cell_types are present in reference_data
  if (!all(cell_types %in% unique(reference_data[[ref_cell_type_col]]))) {
    stop("One or more cell_types specified are not present in reference_data.")
  }
  
  # Function to compute correlation between two cell types
  .computeCorrelation <- function(type1, type2) {
    query_subset <- query_data[ ,query_data[[query_cell_type_col]] == type1, drop = FALSE]
    ref_subset <- reference_data[ ,reference_data[[ref_cell_type_col]] == type2, drop = FALSE]
    
    query_mat <- as.matrix(assay(query_subset, "logcounts"))
    ref_mat <- as.matrix(assay(ref_subset, "logcounts"))
    
    cor_matrix <- cor(query_mat, ref_mat, method = correlation_method)
    mean(cor_matrix)
  }
  
  # Use outer to compute pairwise correlations
  cor_matrix_avg <- outer(cell_types, cell_types, Vectorize(.computeCorrelation))
  
  # Assign cell type names to rows and columns
  dimnames(cor_matrix_avg) <- list(cell_types, cell_types)
  
  return(cor_matrix_avg)
}
