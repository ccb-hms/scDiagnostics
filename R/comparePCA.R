#' @title Compare Principal Components Analysis (PCA) Results
#' 
#' @description This function compares the principal components (PCs) obtained from separate PCA on reference and query 
#' datasets for a single cell type using either cosine similarity or correlation.
#' 
#' @details
#' The function compares the PCs obtained from separate PCA on the reference and query datasets for a single cell type.
#' It computes the similarity between the principal components using either cosine similarity or correlation. 
#' Cosine similarity measures the cosine of the angle between the PC vectors, while correlation computes the 
#' correlation coefficient between them. The function returns a similarity matrix where each element (i, j) 
#' represents the similarity between the i-th PC of the reference dataset and the j-th PC of the query dataset.
#'
#' If cosine similarity is chosen, the function computes the cosine similarity between each pair of PC vectors. 
#' If correlation is chosen, it calculates the correlation coefficient using either Spearman or Pearson method.
#' 
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param n_components Number of principal components to consider. Default is 5.
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
#' @examples
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(ComplexHeatmap)
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
#' # Log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#'
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- getTopHVGs(ref_data, n = 500)
#' query_var <- getTopHVGs(query_data, n = 500)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#'
#' # Subset reference and query data for a specific cell type
#' ref_data_subset <- ref_data_subset[, which(ref_data_subset$reclustered.broad == "CD8")]
#' query_data_subset <- query_data_subset[, which(colData(query_data_subset)$labels == "CD8")]
#'
#' # Run PCA on the reference and query datasets
#' ref_data_subset <- runPCA(ref_data_subset)
#' query_data_subset <- runPCA(query_data_subset)
#'
#' # Call the PCA comparison function
#' similarity_mat <- comparePCA(query_data_subset, ref_data_subset, 
#'                              n_components = 5, 
#'                              metric = c("cosine", "correlation")[1], 
#'                              correlation_method = c("spearman", "pearson")[1])
#'
#' # Create the heatmap
#' heatmap <- Heatmap(similarity_mat,
#'                    name = "Cosine Similarity",
#'                    show_row_names = TRUE,
#'                    show_column_names = TRUE,
#'                    cluster_rows = TRUE,
#'                    cluster_columns = TRUE,
#'                    row_title = paste0("Reference Data (CD4)"), 
#'                    column_title = "Query Data (CD4)", 
#'                    row_dend_side = "left")
#' heatmap
#'
# Compare PCA vectors of reference and query datasets for specific cell type.
comparePCA <- function(reference_data, query_data, 
                       n_components = 5,
                       metric = c("cosine", "correlation")[1], 
                       correlation_method = c("spearman", "pearson")[1]){
    
    # Check if query_data is a SingleCellExperiment object
    if (!is(query_data, "SingleCellExperiment")) {
        stop("query_data must be a SingleCellExperiment object.")
    }
    
    # Check if reference_data is a SingleCellExperiment object
    if (!is(reference_data, "SingleCellExperiment")) {
        stop("reference_data must be a SingleCellExperiment object.")
    }
    
    # Check of genes in both datasets are the same
    if(!all(rownames(attributes(reducedDim(query_data, "PCA"))$rotation) %in%
            rownames(attributes(reducedDim(reference_data, "PCA"))$rotation)))
        stop("The genes in the rotation matrices differ. Consider decreasing the number of genes used for PCA.")
    
    # Check if n_components is a positive integer
    if (!inherits(n_components, "numeric")) {
        stop("n_components should be numeric")
    } else if (any(!n_components == floor(n_components), n_components < 1)) {
        stop("n_components should be an integer, greater than zero.")
    }
    
    # Check if requested number of components is within available components
    if (ncol(reducedDim(reference_data, "PCA")) < n_components) {
        stop("\'n_components\' is larger than number of available components in reference PCA.")
    }
    
    # Check input for metric
    if(!(metric %in% c("cosine", "correlation")))
        stop("\'metric\' should be one of \'cosine\' or \'correlation\'.")
    
    # Check input for correlation method
    if(!(correlation_method %in% c("spearman", "pearson")))
        stop("\'correlation_method\' should be one of \'spearman\' or \'pearson\'.")
    
    # Extract PCA data from reference and query data
    ref_rotation <- attributes(reducedDim(reference_data, "PCA"))$rotation[, 1:n_components]
    query_rotation <- attributes(reducedDim(query_data, "PCA"))$rotation[, 1:n_components]
    
    # Initialize a matrix to store cosine similarities
    similarity_matrix <- matrix(NA, nrow = n_components, ncol = n_components)
    
    if(metric == "cosine"){
        # Function to compute cosine similarity
        .cosine_similarity <- function(x, y) {
            sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
        }
        
        # Loop over each pair of columns and compute cosine similarity
        for (i in 1:n_components) {
            for (j in 1:n_components) {
                similarity_matrix[i, j] <- .cosine_similarity(ref_rotation[, i], query_rotation[, j])
            }
        }
    } else if(metric == "correlation"){
        # Loop over each pair of columns and compute cosine similarity
        for (i in 1:n_components) {
            for (j in 1:n_components) {
                similarity_matrix[i, j] <- cor(ref_rotation[, i], query_rotation[, j], method = correlation_method)
            }
        }
    }
    
    # Add rownames and colnames with % of variance explained for each PC of each dataset 
    rownames(similarity_matrix) <- paste0("Ref PC", 1:n_components, " (", 
                                          round(attributes(reducedDim(reference_data, "PCA"))$varExplained[1:n_components] / 
                                                    sum(attributes(reducedDim(reference_data, "PCA"))$varExplained[1:n_components]) *
                                                            100, 1), "%)")
    colnames(similarity_matrix) <- paste0("Query PC", 1:n_components, " (", 
                                          round(attributes(reducedDim(query_data, "PCA"))$varExplained[1:n_components] / 
                                                    sum(attributes(reducedDim(query_data, "PCA"))$varExplained[1:n_components]) *
                                                    100, 1), "%)")
    
    # Return similarity matrix
    return(similarity_matrix)
}


