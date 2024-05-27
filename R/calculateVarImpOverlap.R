#' @title Compare Gene Importance Across Datasets Using Random Forest
#'
#' @description This function identifies and compares the most important genes for differentiating cell types between a query dataset 
#' and a reference dataset using Random Forest.
#'
#' @details This function uses the Random Forest algorithm to calculate the importance of genes in differentiating between cell types 
#' within both a reference dataset and a query dataset. The function then compares the top genes identified in both datasets to determine 
#' the overlap in their importance scores.
#'
#' @param query_data A SingleCellExperiment object containing the query dataset with logcounts assay.
#' @param reference_data A SingleCellExperiment object containing the reference dataset with logcounts assay.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param n_tree An integer specifying the number of trees to grow in the Random Forest. Default is 500.
#' @param compute_importance A logical value indicating whether to compute variable importance scores for each pairwise combination of 
#' cell types. Default is TRUE.
#' @param n_top An integer specifying the number of top genes to consider when comparing variable importance scores. Default is 20.
#'
#' @return A list containing three elements:
#' \item{var_imp_ref}{A list of data frames containing variable importance scores for each combination of cell types in the reference 
#' dataset.}
#' \item{var_imp_query}{A list of data frames containing variable importance scores for each combination of cell types in the query 
#' dataset.}
#' \item{var_imp_comparison}{A named vector indicating the proportion of top genes that overlap between the reference and query datasets 
#' for each combination of cell types.}
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @examples
#' # Load necessary library
#' library(SingleCellExperiment)
#' library(scuttle)
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
#' # Compare PCA subspaces
#' rf_output <- calculateVarImpOverlap(reference_data = ref_data_subset, query_data = query_data_subset, 
#'                                     query_cell_type_col = "labels", 
#'                                     ref_cell_type_col = "reclustered.broad", 
#'                                     n_tree = 500,
#'                                     n_top = 20)
#' 
# RF function to compare (between datasets) which genes are best at differentiating cell types from each 
calculateVarImpOverlap <- function(query_data, 
                                   reference_data, 
                                   query_cell_type_col, 
                                   ref_cell_type_col,
                                   n_tree = 500,
                                   n_top = 20){
    
    # Extract assay data for reference and query datasets
    ref_x <- t(as.matrix(assay(reference_data, "logcounts")))
    query_x <- t(as.matrix(assay(query_data, "logcounts")))
        
    # Extract labels from reference and query datasets
    ref_y <- reference_data[[ref_cell_type_col]]
    query_y <- query_data[[query_cell_type_col]]
    
    # Remove NA from reference
    ref_x <- ref_x[-which(is.na(ref_y)),]
    ref_y <- ref_y[-which(is.na(ref_y))]
    
    # Finding importance scores for each cell type in reference dataset
    var_imp_ref <- list()
    cell_types <- unique(intersect(ref_y, query_y))
    cell_types_combn <- combn(length(cell_types), 2)
    for(combn_id in 1:ncol(cell_types_combn)){
        
        ref_x_subset <- ref_x[which(ref_y %in% c(cell_types[cell_types_combn[1, combn_id]], cell_types[cell_types_combn[2, combn_id]])),]
        ref_y_subset <- ref_y[which(ref_y %in% c(cell_types[cell_types_combn[1, combn_id]], cell_types[cell_types_combn[2, combn_id]]))]
        training_data <- data.frame(ref_x_subset, cell_type = factor(ref_y_subset))
        rf_binary <- ranger::ranger(cell_type ~ ., data = training_data, num.trees = n_tree, importance = "impurity")
        var_importance_name <- paste0(cell_types[cell_types_combn[1, combn_id]], "-", cell_types[cell_types_combn[2, combn_id]])
        var_imp_ref[[var_importance_name]] <- rf_binary$variable.importance
        var_imp_ref[[var_importance_name]] <- 
            data.frame(Gene = names(var_imp_ref[[var_importance_name]])[order(var_imp_ref[[var_importance_name]], 
                                                                                    decreasing = TRUE)], 
                       RF_Importance = var_imp_ref[[var_importance_name]][order(var_imp_ref[[var_importance_name]], 
                                                                                   decreasing = TRUE)])
    }
    
    # Finding importance scores for each cell type in query dataset
    var_imp_query <- list()
    for(combn_id in 1:ncol(cell_types_combn)){
        
        ref_x_subset <- ref_x[which(ref_y %in% c(cell_types[cell_types_combn[1, combn_id]], cell_types[cell_types_combn[2, combn_id]])),]
        ref_y_subset <- ref_y[which(ref_y %in% c(cell_types[cell_types_combn[1, combn_id]], cell_types[cell_types_combn[2, combn_id]]))]
        training_data <- data.frame(ref_x_subset, cell_type = factor(ref_y_subset))
        rf_binary <- ranger::ranger(cell_type ~ ., data = training_data, num.trees = n_tree, importance = "impurity")
        var_importance_name <- paste0(cell_types[cell_types_combn[1, combn_id]], "-", cell_types[cell_types_combn[2, combn_id]])
        var_imp_query[[var_importance_name]] <- rf_binary$variable.importance
        var_imp_query[[var_importance_name]] <- 
            data.frame(Gene = names(var_imp_query[[var_importance_name]])[order(var_imp_query[[var_importance_name]], 
                                                                                 decreasing = TRUE)], 
                       RF_Importance = var_imp_query[[var_importance_name]][order(var_imp_query[[var_importance_name]], 
                                                                                   decreasing = TRUE)])
    }
    
    # Comparison vector
    var_imp_comparison <- rep(NA, length(var_imp_ref))
    names(var_imp_comparison) <- names(var_imp_ref)
    for(cells in names(var_imp_comparison)){
        var_imp_comparison[cells] <- length(intersect(var_imp_ref[[cells]]$Gene[1:n_top], 
                                                      var_imp_query[[cells]]$Gene[1:n_top])) / n_top
    }
    
    # Return variable importance scores for each combination of cell types in each dataset and the comparison 
    return(list(var_imp_ref = var_imp_ref, 
                var_imp_query = var_imp_query,
                var_imp_comparison = var_imp_comparison))
}