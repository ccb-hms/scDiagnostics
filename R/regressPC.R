#' @title Principal Component Regression
#'
#' @description
#' This function performs linear regression of a covariate of interest onto one
#' or more principal components, based on the data in a \code{\linkS4class{SingleCellExperiment}}
#' object.
#'
#' @details
#' Principal component regression, derived from PCA, can be used to quantify the
#' variance explained by a covariate of interest. Applications for single-cell
#' analysis include quantification of batch effects, assessing clustering
#' homogeneity, and evaluating alignment of query and reference datasets in cell
#' type annotation settings.
#'
#' The function supports multiple regression scenarios:
#' \itemize{
#'   \item Query only, no batch: PC ~ cell_type
#'   \item Query only, with batch: PC ~ cell_type * batch
#'   \item Query + Reference, no batch: PC ~ cell_type * dataset
#'   \item Query + Reference, with batch: PC ~ cell_type * batch (where batch includes Reference)
#' }
#'
#' When batch information is provided with reference data, batches are labeled as
#' "Reference" for reference data and "Query_BatchName" for query batches, with
#' Reference set as the first factor level for interpretation.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' If NULL, the PC scores are regressed against the cell types of the query data.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param query_batch_col The column name in the \code{colData} of \code{query_data} that identifies the batch or sample.
#' If provided, performs interaction analysis with cell types. Default is NULL.
#' @param cell_types A character vector specifying the cell types to include in the analysis. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the analysis. Default is PC1 to PC10.
#' @param adjust_method A character string specifying the method to adjust the p-values.
#'   Options include "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", or "none".
#'   Default is "BH" (Benjamini-Hochberg).
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#' @return
#' A \code{list} containing \itemize{ \item summaries of the linear
#' regression models for each specified principal component, \item the
#' corresponding R-squared (R2) values, \item the variance contributions for
#' each principal component, and \item the total variance explained.}
#'
#' @references Luecken et al. Benchmarking atlas-level data integration in
#' single-cell genomics. Nature Methods, 19:41-50, 2022.
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.regressPCObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Query only analysis
#' regress_res <- regressPC(query_data = query_data,
#'                          query_cell_type_col = "expert_annotation",
#'                          cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                          pc_subset = 1:10)
#' plot(regress_res, plot_type = "r_squared")
#'
#' # Query + Reference analysis
#' regress_res <- regressPC(query_data = query_data,
#'                          reference_data = reference_data,
#'                          query_cell_type_col = "SingleR_annotation",
#'                          ref_cell_type_col = "expert_annotation",
#'                          cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                          pc_subset = 1:10)
#' plot(regress_res, plot_type = "heatmap")
#'
#' @importFrom rlang .data
#' @import SingleCellExperiment
#'
# Function to regress PCs against cell types and batches
regressPC <- function(query_data,
                      reference_data = NULL,
                      query_cell_type_col,
                      ref_cell_type_col = NULL,
                      query_batch_col = NULL,
                      cell_types = NULL,
                      pc_subset = 1:10,
                      adjust_method = c("BH", "holm",
                                        "hochberg", "hommel",
                                        "bonferroni", "BY",
                                        "fdr", "none"),
                      assay_name = "logcounts",
                      max_cells = 2500) {

    # Match argument for adjustment method
    adjust_method <- match.arg(adjust_method)

    # Check standard input arguments
    argumentCheck(reference_data = reference_data,
                  query_data = query_data,
                  ref_cell_type_col = ref_cell_type_col,
                  query_cell_type_col = query_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_query = pc_subset,
                  assay_name = assay_name)

    # Additional check for batch column
    if(!is.null(query_batch_col)){
        if(!query_batch_col %in% colnames(colData(query_data))){
            stop("query_batch_col '", query_batch_col,
                 "' not found in query_data colData.")
        }
    }

    # Downsample reference and query data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    if(!is.null(reference_data)){
        reference_data <- downsampleSCE(sce = reference_data,
                                        max_cells = max_cells)
    }

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        if(is.null(reference_data)){
            cell_types <- na.omit(unique(query_data[[query_cell_type_col]]))
        } else{
            cell_types <- na.omit(unique(c(query_data[[query_cell_type_col]],
                                           reference_data[[ref_cell_type_col]])))
        }
    }

    # Sort cell types for consistent reference category
    cell_types <- sort(cell_types)

    # Perform linear regression for each principal component
    .regressFast <- function(pc, indep_var, df) {
        model <- do.call(speedglm::speedlm,
                         list(formula = paste(pc, " ~ ", indep_var),
                              data = df))
        model_summary <- list(coefficients = summary(model)[["coefficients"]],
                              r_squared = summary(model)[["r.squared"]])
        return(model_summary)
    }

    # Set dependent variables
    dep_vars <- paste0("PC", pc_subset)

    # Case 1: Query data only
    if(is.null(reference_data)){

        # Get query PCA variance for plotting
        query_pca_var <- attr(reducedDim(query_data, "PCA"), "percentVar")

        # Case 1a: No batch - PC ~ cell_type
        if(is.null(query_batch_col)){

            query_labels <- query_data[[query_cell_type_col]]
            regress_data <- data.frame(
                reducedDim(query_data, "PCA")[, dep_vars],
                cell_type = factor(query_data[[query_cell_type_col]],
                                   levels = cell_types))
            regress_data <- regress_data[query_labels %in% cell_types,]

            # Regress PCs against cell types
            summaries <- lapply(dep_vars, .regressFast,
                                indep_var = "cell_type",
                                df = regress_data)
            names(summaries) <- dep_vars

            # Calculate variance contributions
            r_squared <- vapply(summaries, `[[`, numeric(1),
                                x = "r_squared")
            var_expl <- query_pca_var[pc_subset]
            var_contr <- var_expl * r_squared
            total_var_expl <- sum(var_contr)

            regress_res <- list(regression_summaries = summaries,
                                r_squared = r_squared,
                                var_contributions = var_contr,
                                total_variance_explained = total_var_expl,
                                query_pca_var = query_pca_var,
                                indep_var = "cell_type")

            regress_res <- adjustPValues(regress_res,
                                         adjust_method = adjust_method,
                                         indep_var = "cell_type")

            # Case 1b: With batch - PC ~ cell_type * batch
        } else {

            query_labels <- query_data[[query_cell_type_col]]

            regress_data <- data.frame(
                reducedDim(query_data, "PCA")[, dep_vars],
                cell_type = factor(query_data[[query_cell_type_col]],
                                   levels = cell_types),
                batch = factor(query_data[[query_batch_col]]))
            regress_data <- regress_data[query_labels %in% cell_types,]

            # Regress PCs against cell_type * batch interaction
            summaries <- lapply(dep_vars, .regressFast,
                                indep_var = "cell_type * batch",
                                df = regress_data)
            names(summaries) <- dep_vars

            # Calculate variance contributions
            r_squared <- vapply(summaries, `[[`, numeric(1),
                                x = "r_squared")
            var_expl <- query_pca_var[pc_subset]
            var_contr <- var_expl * r_squared
            total_var_expl <- sum(var_contr)

            regress_res <- list(regression_summaries = summaries,
                                r_squared = r_squared,
                                var_contributions = var_contr,
                                total_variance_explained = total_var_expl,
                                query_pca_var = query_pca_var,
                                indep_var = "cell_type_batch_interaction")

            regress_res <- adjustPValues(regress_res,
                                         adjust_method = adjust_method,
                                         indep_var = "cell_type_batch_interaction")
        }

        # Case 2: Query + Reference
    } else {

        # Get reference PCA variance for plotting
        reference_pca_var <- attr(reducedDim(reference_data, "PCA"), "percentVar")

        # Get the projected PCA data
        pca_output <- projectPCA(reference_data = reference_data,
                                 query_data = query_data,
                                 ref_cell_type_col = ref_cell_type_col,
                                 query_cell_type_col = query_cell_type_col,
                                 pc_subset = pc_subset,
                                 assay_name = assay_name)
        pca_output <- pca_output[pca_output[["cell_type"]] %in% cell_types,]

        # Case 2a: No batch - PC ~ cell_type * dataset
        if(is.null(query_batch_col)){

            pca_output[["dataset"]] <- factor(pca_output[["dataset"]],
                                              levels = c("Reference", "Query"))
            pca_output[["cell_type"]] <- factor(pca_output[["cell_type"]],
                                                levels = cell_types)

            # Unified interaction model: PC ~ cell_type * dataset
            summaries <- lapply(dep_vars, .regressFast,
                                indep_var = "cell_type * dataset",
                                df = pca_output)
            names(summaries) <- dep_vars

            # Calculate R-squared
            r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")

            regress_res <- list(regression_summaries = summaries,
                                r_squared = r_squared,
                                reference_pca_var = reference_pca_var,
                                indep_var = "cell_type_dataset_interaction")

            regress_res <- adjustPValues(regress_res,
                                         adjust_method = adjust_method,
                                         indep_var = "cell_type_dataset_interaction")

            # Case 2b: With batch - PC ~ cell_type * batch (where batch includes Reference)
        } else {

            # Create batch labels: "Reference" for reference, "Query_BatchName" for query
            batch_labels <- rep(NA, nrow(pca_output))
            ref_indices <- pca_output[["dataset"]] == "Reference"
            query_indices <- pca_output[["dataset"]] == "Query"

            # Reference gets "Reference" label
            batch_labels[ref_indices] <- "Reference"

            # Query gets "Query_BatchName" labels
            query_cell_names <- rownames(pca_output)[query_indices]
            query_batch_data <- query_data[[query_batch_col]]
            names(query_batch_data) <- colnames(query_data)
            query_batch_values <- query_batch_data[query_cell_names]
            batch_labels[query_indices] <- paste0("Query_", query_batch_values)

            # Add batch information to pca_output
            pca_output[["batch"]] <- factor(
                batch_labels,
                levels = c("Reference",
                           sort(unique(batch_labels[query_indices]))))
            pca_output[["cell_type"]] <- factor(pca_output[["cell_type"]],
                                                levels = cell_types)

            # Remove any rows with missing batch info
            pca_output <- pca_output[!is.na(pca_output[["batch"]]),]

            # Check if this reduces to Case 2a (only one query batch)
            unique_batches <- unique(pca_output[["batch"]])
            if(length(unique_batches) == 2 && "Reference" %in% unique_batches){
                # Fallback to dataset interaction analysis
                pca_output[["dataset"]] <- ifelse(pca_output[["batch"]] == "Reference",
                                                  "Reference", "Query")
                pca_output[["dataset"]] <- factor(pca_output[["dataset"]],
                                                  levels = c("Reference", "Query"))

                # Unified interaction model: PC ~ cell_type * dataset
                summaries <- lapply(dep_vars, .regressFast,
                                    indep_var = "cell_type * dataset",
                                    df = pca_output)
                names(summaries) <- dep_vars

                r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")

                regress_res <- list(regression_summaries = summaries,
                                    r_squared = r_squared,
                                    reference_pca_var = reference_pca_var,
                                    indep_var = "cell_type_dataset_interaction")

                regress_res <- adjustPValues(regress_res,
                                             adjust_method = adjust_method,
                                             indep_var = "cell_type_dataset_interaction")
            } else {
                # True multi-batch case: PC ~ cell_type * batch
                summaries <- lapply(dep_vars, .regressFast,
                                    indep_var = "cell_type * batch",
                                    df = pca_output)
                names(summaries) <- dep_vars

                r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")

                regress_res <- list(regression_summaries = summaries,
                                    r_squared = r_squared,
                                    reference_pca_var = reference_pca_var,
                                    indep_var = "cell_type_batch_interaction")

                regress_res <- adjustPValues(regress_res,
                                             adjust_method = adjust_method,
                                             indep_var = "cell_type_batch_interaction")
            }
        }
    }

    # Return regression output
    class(regress_res) <- c(class(regress_res), "regressPCObject")
    return(regress_res)
}

#' @title Adjust P-Values in Regression Results
#'
#' @description
#' Adjusts the p-values in the regression results using a specified adjustment method.
#' The adjustment is performed for different regression types including cell type,
#' dataset, and cell type-batch interaction analyses.
#'
#' @details
#' This function adjusts p-values from regression results stored in a list. The adjustment
#' can be applied across different regression structures depending on the analysis type.
#' The method for adjusting p-values can be selected from various options such as
#' Benjamini-Hochberg (BH), Holm, and others, which are supported by the `p.adjust` function in R.
#'
#' @param regress_res A list containing regression results. The structure varies by analysis type.
#' @param adjust_method A character string specifying the method to adjust the p-values.
#'   Options include "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", or "none".
#'   Default is "BH" (Benjamini-Hochberg).
#' @param indep_var A character string specifying the independent variable for the adjustment.
#'   Options are "cell_type", "cell_type_dataset_interaction", or "cell_type_batch_interaction".
#'
#' @keywords internal
#'
#' @return A list similar to \code{regress_res}, but with an added column for adjusted p-values
#'   in the coefficients tables.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @importFrom stats p.adjust
#'
# Function to compute adjusted p-values for PC regression
adjustPValues <- function(regress_res,
                          adjust_method = c("BH", "holm",
                                            "hochberg", "hommel",
                                            "bonferroni", "BY",
                                            "fdr", "none"),
                          indep_var = c("cell_type", "cell_type_dataset_interaction",
                                        "cell_type_batch_interaction")){

    # Match arguments
    adjust_method <- match.arg(adjust_method)
    indep_var <- match.arg(indep_var)

    # All interaction models have the same structure
    if(indep_var %in% c("cell_type", "cell_type_dataset_interaction", "cell_type_batch_interaction")){
        # Add adjusted p-values for unified model structure
        for (pc in names(regress_res[["regression_summaries"]])) {
            coeffs <- regress_res[["regression_summaries"]][[pc]][["coefficients"]]
            p.adjusted <- p.adjust(coeffs[["p.value"]], method = adjust_method)
            coeffs[["p.adjusted"]] <- p.adjusted
            regress_res[["regression_summaries"]][[pc]][["coefficients"]] <- coeffs
        }
    }

    # Return regression object
    return(regress_res)
}
