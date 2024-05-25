#' @title Random Forest Cell Annotation and Variable Importance Scores
#'
#' @description Perform Random Forest Diagnostics on Query and Reference Single-Cell RNA-Seq Data.
#'
#' @details This function trains a Random Forest (RF) model using a reference dataset and then predicts cell types in a query dataset. 
#' It provides a comparison table of true vs. predicted cell types. Optionally, it can also compute variable importance scores for each 
#' possible pairwise combination of cell types, providing insights into the genes that are most important for distinguishing between 
#' each pair of cell types.
#'
#' @param query_data A SingleCellExperiment object containing the query dataset with logcounts assay.
#' @param reference_data A SingleCellExperiment object containing the reference dataset with logcounts assay.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param n_tree An integer specifying the number of trees to grow in the Random Forest. Default is 500.
#' @param compute_importance A logical value indicating whether to compute variable importance scores for each pairwise combination of cell types. Default is TRUE.
#'
#'
#' @return A list with the following components:
#' \describe{
#'   \item{rf_pred}{A character vector of predicted cell types for the query dataset.}
#'   \item{comparison_table}{A table comparing true cell annotations and RF predictions for the query dataset.}
#'   \item{var_importance}{A list of data frames containing the variable importance scores for each pair of cell types. Only returned if \code{compute_importance} is TRUE.}
#' }
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
#' sce <- scRNAseq::HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
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
#' ref_var <- scran::getTopHVGs(ref_data, n = 500)
#' query_var <- scran::getTopHVGs(query_data, n = 500)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#' 
#' # Compare PCA subspaces
#' rf_output <- RFDiagnostics(reference_data = ref_data_subset, query_data = query_data_subset, 
#'                            query_cell_type_col = "labels", 
#'                            ref_cell_type_col = "reclustered.broad", 
#'                            n_tree = 500)
#' 
# RF function to get cell annotations and compute variable importance scores
RFDiagnostics <- function(query_data, 
                          reference_data, 
                          query_cell_type_col, 
                          ref_cell_type_col,
                          n_tree = 500,
                          compute_importance = TRUE){
    
    # Extract assay data for reference and query datasets
    ref_x <- t(as.matrix(assay(reference_data, "logcounts")))
    query_x <- t(as.matrix(assay(query_data, "logcounts")))
        
    # Extract labels from reference and query datasets
    ref_y <- reference_data[[ref_cell_type_col]]
    query_y <- query_data[[query_cell_type_col]]
    
    # Remove NA from reference
    ref_x <- ref_x[-which(is.na(ref_y)),]
    ref_y <- ref_y[-which(is.na(ref_y))]
    
    # Train the RF model 
    training_data <- data.frame(ref_x, cell_type = factor(ref_y))
    rf_multiclass <- ranger::ranger(cell_type ~ ., data = training_data, num.trees = n_tree)

    # Predict the class
    rf_pred <- predict(rf_multiclass, data.frame(query_x))
    comparison_table <- table(Cell_Annotation = query_y, RF_Prediction = rf_pred$predictions)
    
    # Finding importance scores for each cell type
    var_importance <- list()
    cell_types <- unique(query_y)
    cell_types_combn <- combn(length(cell_types), 2)
    for(combn_id in 1:ncol(cell_types_combn)){
        
        ref_x_subset <- ref_x[which(ref_y %in% c(cell_types[cell_types_combn[1, combn_id]], cell_types[cell_types_combn[2, combn_id]])),]
        ref_y_subset <- ref_y[which(ref_y %in% c(cell_types[cell_types_combn[1, combn_id]], cell_types[cell_types_combn[2, combn_id]]))]
        training_data <- data.frame(ref_x_subset, cell_type = factor(ref_y_subset))
        rf_binary <- ranger::ranger(cell_type ~ ., data = training_data, num.trees = n_tree, importance = "impurity")
        var_importance_name <- paste0(cell_types[cell_types_combn[1, combn_id]], "-", cell_types[cell_types_combn[2, combn_id]])
        var_importance[[var_importance_name]] <- rf_binary$variable.importance
        var_importance[[var_importance_name]] <- 
            data.frame(Gene = names(var_importance[[var_importance_name]])[order(var_importance[[var_importance_name]], 
                                                                                    decreasing = TRUE)], 
                       RF_Importance = var_importance[[var_importance_name]][order(var_importance[[var_importance_name]], 
                                                                                   decreasing = TRUE)])
    }
    
    # Return RF predictions, the prediction comparison table and variable importance scores
    return(list(rf_pred = as.character(rf_pred$predictions), 
                comparison_table = comparison_table,
                var_importance = var_importance))
}