#' @title Compare Subspaces Spanned by Top Principal Components Using Canonical Correlation Analysis
#' 
#' @description
#' This function compares the subspaces spanned by the top principal components (PCs) of the reference 
#' and query datasets using canonical correlation analysis (CCA). It calculates the canonical variables, 
#' correlations, and a similarity measure for the subspaces.
#'
#' @details
#' This function performs canonical correlation analysis (CCA) to compare the subspaces spanned by the 
#' top principal components (PCs) of the reference and query datasets. The function extracts the rotation 
#' matrices corresponding to the specified PCs and performs CCA on these matrices. It computes the canonical 
#' variables and their corresponding correlations. Additionally, it calculates a similarity measure for the 
#' canonical variables using cosine similarity. The output is a list containing the canonical coefficients 
#' for both datasets, the cosine similarity values, and the canonical correlations.

#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) 
#' to compare. Default is the first five PCs.
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
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.compareCCA}}
#' 
#' @examples
#' \donttest{
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(ggplot2)
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
#' ref_data_subset <- runPCA(ref_data_subset, ncomponents = 50)
#' query_data_subset <- runPCA(query_data_subset, ncomponents = 50)
#' 
#' # Compare CCA
#' cca_comparison <- compareCCA(query_data_subset, ref_data_subset, 
#'                              pc_subset = c(1:5))
#' 
#' # Visualize output of CCA comparison
#' plot(cca_comparison)
#' }
#' 
# Function to compare subspace spanned by top PCs in reference and query datasets
compareCCA <- function(reference_data, query_data, 
                       pc_subset = c(1:5)){
    
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
        stop("The genes in the rotation matrices differ. Consider decreasing the number of genes using for PCA.")
    
    # Check input if PC subset is valid
    if(!all(c(pc_subset %in% 1:ncol(reducedDim(reference_data, "PCA")), 
              pc_subset %in% 1:ncol(reducedDim(query_data, "PCA")))))
        stop("\'pc_subset\' is out of range.")
    
    # Extract the rotation matrices
    ref_pcs <- attributes(reducedDim(reference_data, "PCA"))$rotation[, pc_subset]
    query_pcs <- attributes(reducedDim(query_data, "PCA"))$rotation[, pc_subset]
    
    # Perform CCA
    cca_result <- cancor(ref_pcs, query_pcs)
    
    # Extract canonical variables and correlations
    canonical_ref <- cca_result$xcoef
    canonical_query <- cca_result$ycoef
    correlations <- cca_result$cor
    
    # Function to compute similarity measure (e.g., cosine similarity)
    .cosine_similarity <- function(u, v) {
        return(abs(sum(u * v)) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
    }
    
    # Compute similarities and account for correlations
    similarities <- rep(0, length(pc_subset))
    for (i in 1:length(pc_subset)) {
        similarities[i] <- .cosine_similarity(canonical_ref[, i], canonical_query[, i])
    }
    
    # Update class of return output
    output <- list(coef_ref = canonical_ref,
                   coef_query = canonical_query,
                   cosine_similarity = similarities,
                   correlations = correlations)
    class(output) <- c(class(output), "compareCCA")

    # Return cosine similarity output
    return(output)
}

