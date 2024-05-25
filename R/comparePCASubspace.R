#' @title Compare Subspaces Spanned by Top Principal Components
#' 
#' @description
#' This function compares the subspaces spanned by the top principal components of 
#' two datasets using a weighted cosine similarity measure. It computes the cosine 
#' similarity between the principal angles of the subspaces, weighted by the proportion 
#' of variance explained by each principal component.
#'
#' @details
#' The function takes as input two SingleCellExperiment objects representing the 
#' reference and query datasets, respectively. Users can specify a subset of principal 
#' components (PCs) to compare. The cosine similarity between the principal angles of 
#' the subspaces spanned by the specified PCs is computed and weighted by the proportion 
#' of variance explained by each PC.
#'
#' @param reference_data A SingleCellExperiment object representing the reference dataset.
#' @param query_data A SingleCellExperiment object representing the query dataset.
#' @param pc_subset A numeric vector specifying the subset of principal components (PCs) 
#' to compare. Default is the first five PCs.
#'
#' @return A list containing the following components:
#'   \item{principal_angles_cosines}{A numeric vector of cosine values of principal angles.}
#'   \item{average_variance_explained}{A numeric vector of average variance explained by each PC.}
#'   \item{weighted_cosine_similarity}{A numeric value representing the weighted cosine similarity.}
#'
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @examples
#' # Load necessary library
#' library(SingleCellExperiment)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(ggplot2)
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
#' # Compare PCA subspaces
#' subspace_comparison <- comparePCASubspace(reference_data = ref_data_subset, 
#'                                           query_data = query_data_subset, 
#'                                           pc_subset = c(1:5))
#' 
#' # Create a data frame for plotting
#' comparison_data <- data.frame(PC = paste0("PC", pc_subset),
#'                               Cosine = subspace_comparison$principal_angles_cosines,
#'                               VarianceExplained = subspace_comparison$average_variance_explained)
#' 
#' # Plot the cosines of principal angles with variance explained as size
#' ggplot(data, aes(x = PC, y = Cosine, size = VarianceExplained)) +
#'     geom_point() +
#'     scale_size_continuous(range = c(3, 10)) +
#'     labs(title = "Principal Angles Cosines with Variance Explained",
#'          x = "Principal Component",
#'          y = "Cosine of Principal Angle",
#'          size = "Variance Explained") +
#'     theme_minimal()
#' 
# Function to compare subspace spanned by top PCs in reference and query datasets
comparePCASubspace <- function(reference_data, query_data, 
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
    
    # Compute the inner product of the selected top principal components
    inner_product <- t(ref_pcs) %*% query_pcs
    
    # Compute the cosine similarity (cosine of principal angle)
    cosine_similarity <- abs(inner_product) / (sqrt(apply(ref_pcs^2, 2, sum)) * sqrt(apply(query_pcs^2, 2, sum)))
    
    # Vector to store top cosine similarities
    top_cosine <- numeric(length(pc_subset))
    # Matrix to store PC IDs for each top cosine similarity
    cosine_id <- matrix(NA, nrow = length(pc_subset), ncol = 2)
    colnames(cosine_id) <- c("Ref", "Query")
    
    # Looping to store top cosine similarities and PC IDs
    for(id in 1:length(pc_subset)){
        
        # Store data for top cosine
        top_ref <- which.max(apply(cosine_similarity, 1, max))
        top_query <- which.max(cosine_similarity[top_ref,])
        top_cosine[id] <- cosine_similarity[top_ref, top_query]
        cosine_id[id,] <- c(top_ref, top_query)
        
        # Remove as candidate
        cosine_similarity[top_ref,] <- -Inf 
        cosine_similarity[, top_query] <- -Inf
    }
    
    # Vector of variance explained
    var_explained_ref <- attributes(reducedDim(reference_data, "PCA"))$varExplained[pc_subset]/
        sum(attributes(reducedDim(reference_data, "PCA"))$varExplained[pc_subset])
    var_explained_query <- attributes(reducedDim(reference_data, "PCA"))$varExplained[pc_subset]/
        sum(attributes(reducedDim(reference_data, "PCA"))$varExplained[pc_subset])
    var_explained_avg <- (var_explained_ref[cosine_id[, 1]] + var_explained_query[cosine_id[, 2]]) / 2
    
    # Weighted cosine similarity score
    weighted_cosine_similarity <- sum(top_cosine * var_explained_avg)
    
    # Return cosine similarity output
    return(list(cosine_similarity = top_cosine,
                cosine_id = cosine_id,
                var_explained_avg = var_explained_avg,
                weighted_cosine_similarity = weighted_cosine_similarity))
}
