#' @title Compare Gene Importance Across Datasets Using Random Forest
#'
#' @description
#' This function identifies and compares the most important genes for differentiating cell types between a query dataset
#' and a reference dataset using Random Forest.
#'
#' @details This function uses the Random Forest algorithm to calculate the importance of genes in differentiating between cell types
#' within both a reference dataset and a query dataset. The function then compares the top genes identified in both datasets to determine
#' the overlap in their importance scores.
#'
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If NULL, then the variable importance scores are only computed for the reference data. Default is NULL.
#' @param ref_cell_type_col A character string specifying the column name in the reference dataset containing cell type annotations.
#' @param query_cell_type_col A character string specifying the column name in the query dataset containing cell type annotations.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param n_tree An integer specifying the number of trees to grow in the Random Forest. Default is 500.
#' @param n_top An integer specifying the number of top genes to consider when comparing variable importance scores. Default is 50.
#' @param assay_name Name of the assay on which to perform computations. Defaults to \code{"logcounts"}.
#'
#' @return A list containing three elements:
#' \item{var_imp_ref}{A list of data frames containing variable importance scores for each combination of cell types in the reference
#' dataset.}
#' \item{var_imp_query}{A list of data frames containing variable importance scores for each combination of cell types in the query
#' dataset.}
#' \item{var_imp_comparison}{A named vector indicating the proportion of top genes that overlap between the reference and query
#' datasets for each combination of cell types.}
#'
#' @export
#'
#' @references
#' Breiman, L. (2001). "Random forests". *Machine Learning*, 45(1), 5-32. doi:10.1023/A:1010933404324.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compute important variables for all pairwise cell comparisons
#' rf_output <- calculateVarImpOverlap(reference_data = reference_data,
#'                                     query_data = query_data,
#'                                     query_cell_type_col = "SingleR_annotation",
#'                                     ref_cell_type_col = "expert_annotation",
#'                                     n_tree = 500,
#'                                     n_top = 50)
#'
#' # Comparison table
#' rf_output$var_imp_comparison
#'
# RF function to compare (between datasets) which genes are best at differentiating cell types from each
calculateVarImpOverlap <- function(reference_data,
                                   query_data = NULL,
                                   ref_cell_type_col,
                                   query_cell_type_col = NULL,
                                   cell_types = NULL,
                                   n_tree = 500,
                                   n_top = 50,
                                   assay_name = "logcounts"){

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types)

    # Check if n_tree is a positive integer
    if (!is.numeric(n_tree) || n_tree <= 0 || n_tree != as.integer(n_tree)) {
        stop("\'n_tree\' must be a positive integer.")
    }

    # Check if n_top is a positive integer
    if (!is.numeric(n_top) || n_top <= 0 || n_top != as.integer(n_top)) {
        stop("\'n_top\' must be a positive integer.")
    }

    # Get cell types if they are not specified by user
    if(is.null(cell_types)){
        if(is.null(query_data)){
            cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]])))
        } else{
            cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                           query_data[[query_cell_type_col]])))
        }
    }

    # Extract assay data for reference and query datasets
    ref_x <- t(as.matrix(assay(reference_data, assay_name)))

    # Extract labels from reference and query datasets
    ref_y <- reference_data[[ref_cell_type_col]]

    # Remove NA from reference
    ref_x <- ref_x[which(!is.na(ref_y)),]
    ref_y <- na.omit(ref_y)

    # Finding importance scores for each cell type in reference dataset
    var_imp_ref <- list()
    cell_types <- unique(ref_y)
    cell_types_combn <- combn(length(cell_types), 2)
    for(combn_id in seq_len(ncol(cell_types_combn))){

        ref_x_subset <- ref_x[which(
            ref_y %in% c(cell_types[cell_types_combn[1, combn_id]],
                         cell_types[cell_types_combn[2, combn_id]])),]
        ref_y_subset <- ref_y[which(
            ref_y %in% c(cell_types[cell_types_combn[1, combn_id]],
                         cell_types[cell_types_combn[2, combn_id]]))]
        training_data <- data.frame(ref_x_subset,
                                    cell_type = factor(ref_y_subset))
        rf_binary <- ranger::ranger(cell_type ~ ., data = training_data,
                                    num.trees = n_tree, importance = "impurity")
        var_importance_name <- paste0(cell_types[cell_types_combn[1, combn_id]],
                                      "-",
                                      cell_types[cell_types_combn[2, combn_id]])
        var_imp_ref[[var_importance_name]] <- rf_binary$variable.importance
        names(var_imp_ref[[var_importance_name]]) <- colnames(ref_x_subset)
        var_imp_ref[[var_importance_name]] <-
            data.frame(Gene = names(
                var_imp_ref[[var_importance_name]])[order(
                    var_imp_ref[[var_importance_name]],
                    decreasing = TRUE)],
                       RF_Importance = var_imp_ref[[var_importance_name]][order(
                           var_imp_ref[[var_importance_name]],
                           decreasing = TRUE)])
        rownames(var_imp_ref[[var_importance_name]]) <- NULL
    }

    if(!is.null(query_data)){

        # Extract assay data for query dataset
        query_x <- t(as.matrix(assay(query_data, assay_name)))

        # Extract labels from query dataset
        query_y <- query_data[[query_cell_type_col]]

        # Remove NA from query
        query_x <- query_x[which(!is.na(query_y)),]
        query_y <- na.omit(query_y)

        # Finding importance scores for each cell type in query dataset
        cell_types <- unique(intersect(ref_y, query_y))
        var_imp_query <- list()
        for(combn_id in seq_len(ncol(cell_types_combn))){

            query_x_subset <- query_x[which(
                query_y %in% c(cell_types[cell_types_combn[1, combn_id]],
                               cell_types[cell_types_combn[2, combn_id]])),]
            query_y_subset <- query_y[which(
                query_y %in% c(cell_types[cell_types_combn[1, combn_id]],
                               cell_types[cell_types_combn[2, combn_id]]))]
            training_data <- data.frame(query_x_subset,
                                        cell_type = factor(query_y_subset))
            rf_binary <- ranger::ranger(cell_type ~ ., data = training_data,
                                        num.trees = n_tree,
                                        importance = "impurity")
            var_importance_name <- paste0(
                cell_types[cell_types_combn[1, combn_id]], "-",
                cell_types[cell_types_combn[2, combn_id]])
            var_imp_query[[var_importance_name]] <-
                rf_binary$variable.importance
            names(var_imp_query[[var_importance_name]]) <-
                colnames(query_x_subset)
            var_imp_query[[var_importance_name]] <-
                data.frame(Gene = names(
                    var_imp_query[[var_importance_name]])[order(
                        var_imp_query[[var_importance_name]],
                        decreasing = TRUE)],
                           RF_Importance = var_imp_query[[var_importance_name]][order(
                               var_imp_query[[var_importance_name]],
                               decreasing = TRUE)])
            rownames(var_imp_query[[var_importance_name]]) <- NULL
        }

        # Comparison vector
        var_imp_comparison <- rep(NA, length(var_imp_ref))
        names(var_imp_comparison) <- names(var_imp_ref)
        for(cells in names(var_imp_comparison)){
            var_imp_comparison[cells] <- length(
                intersect(var_imp_ref[[cells]][["Gene"]][seq_len(n_top)],
                          var_imp_query[[cells]][["Gene"]][seq_len(n_top)])) /
                n_top
        }

        # Return variable importance scores for each combination of cell types in each dataset and the comparison
        return(list(var_imp_ref = var_imp_ref,
                    var_imp_query = var_imp_query,
                    var_imp_comparison = var_imp_comparison))
    } else{
        return(list(var_imp_ref = var_imp_ref))
    }
}

