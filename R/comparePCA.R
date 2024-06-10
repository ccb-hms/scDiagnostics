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
#' @seealso \code{\link{plot.comparePCA}}
#' 
#' @examples
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(scater)
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
#'                              pc_subset = c(1:5), 
#'                              n_top_vars = 50,
#'                              metric = c("cosine", "correlation")[1], 
#'                              correlation_method = c("spearman", "pearson")[1])
#'
#' # Create the heatmap
#' plot(similarity_mat)
#' 
# Compare PCA vectors of reference and query datasets for specific cell type.
comparePCA <- function(reference_data, query_data, 
                       pc_subset = c(1:5),
                       n_top_vars = 50,
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
    
    # Check input if PC subset is valid
    if(!all(c(pc_subset %in% 1:ncol(reducedDim(reference_data, "PCA")), 
              pc_subset %in% 1:ncol(reducedDim(query_data, "PCA")))))
        stop("\'pc_subset\' is out of range.")
    
    # Check input for metric
    if(!(metric %in% c("cosine", "correlation")))
        stop("\'metric\' should be one of \'cosine\' or \'correlation\'.")
    
    # Check input for correlation method
    if(!(correlation_method %in% c("spearman", "pearson")))
        stop("\'correlation_method\' should be one of \'spearman\' or \'pearson\'.")
    
    # Extract PCA data from reference and query data
    ref_rotation <- attributes(reducedDim(reference_data, "PCA"))$rotation[, pc_subset]
    query_rotation <- attributes(reducedDim(query_data, "PCA"))$rotation[, pc_subset]
    
    # Function to identify high-loading variables for each PC
    .getHighLoadingVars <- function(rotation_mat, n_top_vars) {
        high_loading_vars <- lapply(1:ncol(rotation_mat), function(pc) {
            abs_loadings <- abs(rotation_mat[, pc])
            top_vars <- names(sort(abs_loadings, decreasing = TRUE))[1:n_top_vars]
            return(top_vars)
        })
        return(high_loading_vars)
    }
    
    # Get union of variables with highest loadings
    top_ref <- .getHighLoadingVars(ref_rotation, n_top_vars)
    top_query <- .getHighLoadingVars(query_rotation, n_top_vars)
    top_union <- lapply(1:length(pc_subset), function(i) return(union(top_ref[[i]], top_query[[i]])))

    # Initialize a matrix to store cosine similarities
    similarity_matrix <- matrix(NA, nrow = length(pc_subset), ncol = length(pc_subset))
    
    if(metric == "cosine"){
        # Function to compute cosine similarity
        .cosine_similarity <- function(x, y) {
            sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
        }
        
        # Loop over each pair of columns and compute cosine similarity
        for (i in 1:length(pc_subset)) {
            for (j in 1:length(pc_subset)) {
                combination_union <- union(top_union[[i]], top_union[[j]])
                similarity_matrix[i, j] <- .cosine_similarity(ref_rotation[combination_union, i], query_rotation[combination_union, j])
            }
        }
    } else if(metric == "correlation"){
        # Loop over each pair of columns and compute cosine similarity
        for (i in 1:length(pc_subset)) {
            for (j in 1:length(pc_subset)) {
                combination_union <- union(top_union[[i]], top_union[[j]])
                similarity_matrix[i, j] <- cor(ref_rotation[combination_union, i], query_rotation[combination_union, j], 
                                               method = correlation_method)
            }
        }
    }
    
    # Add rownames and colnames with % of variance explained for each PC of each dataset 
    rownames(similarity_matrix) <- paste0("Ref PC", pc_subset, " (", 
                                          round(attributes(reducedDim(reference_data, "PCA"))$varExplained[pc_subset] / 
                                                    sum(attributes(reducedDim(reference_data, "PCA"))$varExplained[pc_subset]) *
                                                            100, 1), "%)")
    colnames(similarity_matrix) <- paste0("Query PC", pc_subset, " (", 
                                          round(attributes(reducedDim(query_data, "PCA"))$varExplained[pc_subset] / 
                                                    sum(attributes(reducedDim(query_data, "PCA"))$varExplained[pc_subset]) *
                                                    100, 1), "%)")
    
    # Update class of return output
    class(similarity_matrix) <- c(class(similarity_matrix), "comparePCA")
    
    # Return similarity matrix
    return(similarity_matrix)
}


