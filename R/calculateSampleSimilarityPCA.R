#' @title Calculate Sample Similarity Using PCA Loadings
#'
#' @description 
#' This function calculates the cosine similarity between samples based on the principal components (PCs)
#' obtained from PCA (Principal Component Analysis) loadings.
#'
#' @details 
#' This function calculates the cosine similarity between samples based on the loadings of the selected
#' principal components obtained from PCA. It extracts the rotation matrix from the PCA results of the 
#' \code{\linkS4class{SingleCellExperiment}} object and identifies the high-loading variables for each selected PC. 
#' Then, it computes the cosine similarity between samples using the high-loading variables for each PC.
#'
#' @param se_object A \code{\linkS4class{SingleCellExperiment}} object containing expression data.
#' @param samples A character vector specifying the samples for which to compute the similarity.
#' @param pc_subset A numeric vector specifying the subset of principal components to consider (default: c(1:5)).
#' @param n_top_vars An integer indicating the number of top loading variables to consider for each PC (default: 50).
#'
#' @return A data frame containing cosine similarity values between samples for each selected principal component.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.calculateSampleSimilarityPCA}}
#'
#' @examples
#' # Load required libraries
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scran)
#' library(scater)
#'
#' # Load data (replace with your data loading)
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#' 
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#' 
#' # log transform datasets
#' ref_data <- scuttle::logNormCounts(ref_data)
#' query_data <- scuttle::logNormCounts(query_data)
#' 
#' # Get cell type scores using SingleR (or any other cell type annotation method)
#' scores <- SingleR::SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#' 
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Selecting highly variable genes (can be customized by the user)
#' ref_var <- scran::getTopHVGs(ref_data, n = 2000)
#' query_var <- scran::getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#'
#' # Run PCA on the reference data (assumed to be prepared)
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' # Store PCA anomaly data and plots
#' anomaly_output <- detectAnomaly(reference_data = ref_data_subset, 
#'                                 ref_cell_type_col = "reclustered.broad", 
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5) 
#' top6_anomalies <- names(sort(anomaly_output$Combined$reference_anomaly_scores, decreasing = TRUE)[1:6])
#' 
#' # Compute cosine similarity between anomalies and top PCs
#' cosine_similarities <- calculateSampleSimilarityPCA(ref_data_subset, samples = top6_anomalies, 
#'                                                     pc_subset = c(1:10), n_top_vars = 50)
#' cosine_similarities
#' 
#' # Plot similarities
#' plot(cosine_similarities, pc_subset = c(1:5))
#' 
# Function to calculate cosine similarities between samples and PCs
calculateSampleSimilarityPCA <- function(se_object, samples, pc_subset = c(1:5), n_top_vars = 50){
    
    # Extract rotation matrix for SingleCellExperiment object
    rotation_mat <- attributes(reducedDim(se_object, "PCA"))$rotation[, pc_subset]
    
    # FUnction to identify high-loading variables for each PC
    .getHighLoadingVars <- function(rotation_mat, n_top_vars) {
        high_loading_vars <- lapply(1:ncol(rotation_mat), function(pc) {
            abs_loadings <- abs(rotation_mat[, pc])
            top_vars <- names(sort(abs_loadings, decreasing = TRUE))[1:n_top_vars]
            return(top_vars)
        })
        return(high_loading_vars)
    }
    
    # Extract high-loading variables
    high_loading_vars <- .getHighLoadingVars(rotation_mat, n_top_vars)
    
    # Function to compute cosine similarity
    .cosine_similarity <- function(vector1, vector2) {
        sum(vector1 * vector2) / (sqrt(sum(vector1^2)) * sqrt(sum(vector2^2)))
    }
    
    # Function to compute cosine similarity for each PC using high-loading variables
    .computeCosineSimilarity <- function(samples, rotation_mat, high_loading_vars) {
        similarities <- lapply(1:length(high_loading_vars), function(pc) {
            vars <- high_loading_vars[[pc]]
            sample_subset <- samples[, vars, drop = FALSE]
            pc_vector <- rotation_mat[vars, pc]
            apply(sample_subset, 1, .cosine_similarity, vector2 = pc_vector)
        })
        return(similarities)
    }
    
    # Calculate similarities
    assay_mat <- t(as.matrix(assay(se_object[, samples], "logcounts")))
    similarities <- .computeCosineSimilarity(assay_mat, rotation_mat, high_loading_vars)
    
    # Format the result into a data frame for easy interpretation
    similarity_df <- do.call(cbind, similarities)
    colnames(similarity_df) <- paste0("PC", 1:ncol(rotation_mat))
    
    # Update class of output
    class(similarity_df) <- c(class(similarity_df), "calculateSampleSimilarityPCA")
    return(similarity_df)
}