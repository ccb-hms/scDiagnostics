#' @title Perform Hotelling's T-squared Test on PCA Scores for Single-cell RNA-seq Data
#'
#' @description This function performs Hotelling's T-squared test to assess the similarity between reference and query datasets 
#' for each cell type based on their PCA scores.
#'
#' @details This function first performs PCA on the reference dataset and then projects the query dataset onto the PCA space 
#' of the reference data. For each cell type, it computes pseudo-bulk signatures in the PCA space by averaging the principal 
#' component scores of cells belonging to that cell type. Hotelling's T-squared test is then performed to compare the mean 
#' vectors of the pseudo-bulk signatures between the reference and query datasets. The resulting p-values indicate the similarity 
#' between the reference and query datasets for each cell type.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param n_components An integer specifying the number of principal components to use for projection. Defaults to 10. 
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#'
#' @return A named numeric vector of p-values from Hotelling's T-squared test for each cell type.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
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
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset, ncomponents = 50)
#'
#' # Get the p-values from the test
#' p_values <- calculateHotellingPValue(query_data_subset, ref_data_subset, 
#'                                      n_components = 10, 
#'                                      query_cell_type_col = "labels", 
#'                                      ref_cell_type_col = "reclustered.broad",
#'                                      pc_subset = c(1:10)) 
#' round(p_values, 5)
#'                          
# Function to perform Hotelling T^2 test for each cell type
# The test is performed on the PCA space of the reference data 
# The query data projected onto PCA space of reference
calculateHotellingPValue <- function(query_data, reference_data, 
                                     n_components = 10, 
                                     query_cell_type_col, 
                                     ref_cell_type_col,
                                     pc_subset = c(1:5)) {
    
    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data, reference_data = reference_data, 
                             n_components = n_components, 
                             query_cell_type_col = query_cell_type_col, 
                             ref_cell_type_col = ref_cell_type_col, 
                             return_value = "list")
    
    # Get unique cell types
    unique_cell_types <- na.omit(unique(c(colData(reference_data)[[ref_cell_type_col]],
                                          colData(query_data)[[query_cell_type_col]])))
    
    # Create a list to store p-values for each cell type
    p_values <- rep(NA, length(unique_cell_types))
    names(p_values) <- unique_cell_types
    
    for (cell_type in unique_cell_types) {
        
        # Subset principal component scores for current cell type
        ref_subset_scores <- pca_output$ref[which(cell_type == reference_data[[ref_cell_type_col]]), pc_subset]
        query_subset_scores <- pca_output$query[which(cell_type == query_data[[query_cell_type_col]]), pc_subset]
        # Calculate the p-value
        hotelling_output <- Hotelling::hotelling.test(x = ref_subset_scores, y = query_subset_scores)
        # Store the result
        p_values[cell_type] <- hotelling_output$pval
    }
    
    # Return p-values
    return(p_values)
}
