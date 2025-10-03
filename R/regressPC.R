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
#'   \item Query only, no batch: PC  cell_type
#'   \item Query only, with batch: PC  cell_type * batch
#'   \item Query + Reference, no batch: PC  cell_type * dataset
#'   \item Query + Reference, with batch: PC  cell_type * batch (where batch includes Reference)
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
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is 5000.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is 5000.
#'
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
#' # Visualize results
#' plot(regress_res, plot_type = "r_squared")
#' plot(regress_res, plot_type = "variance_contribution")
#' plot(regress_res, plot_type = "coefficient_heatmap")
#'
#' # Query + Reference analysis
#' regress_res <- regressPC(query_data = query_data,
#'                          reference_data = reference_data,
#'                          query_cell_type_col = "SingleR_annotation",
#'                          ref_cell_type_col = "expert_annotation",
#'                          cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                          pc_subset = 1:10)
#' # Visualize results
#' plot(regress_res, plot_type = "r_squared")
#' plot(regress_res, plot_type = "variance_contribution")
#' plot(regress_res, plot_type = "coefficient_heatmap")
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
                      max_cells_ref = 5000,
                      max_cells_query = 5000) {

    # Match argument for adjustment method
    adjust_method <- match.arg(adjust_method)

    # Check standard input arguments
    argumentCheck(reference_data = reference_data,
                  query_data = query_data,
                  ref_cell_type_col = ref_cell_type_col,
                  query_cell_type_col = query_cell_type_col,
                  pc_subset_query = pc_subset,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    query_cols_to_convert <- c(query_cell_type_col)
    if (!is.null(query_batch_col)) {
        query_cols_to_convert <- c(query_cols_to_convert, query_batch_col)
    }
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cols_to_convert)

    if(!is.null(reference_data)){
        reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                    convert_cols = ref_cell_type_col)
    }

    # Additional check for batch column
    if(!is.null(query_batch_col)){
        if(!query_batch_col %in% colnames(colData(query_data))){
            stop("query_batch_col '", query_batch_col,
                 "' not found in query_data colData.")
        }
    }

    # Select cell types
    cell_types <- selectCellTypes(query_data = query_data,
                                  reference_data = reference_data,
                                  query_cell_type_col = query_cell_type_col,
                                  ref_cell_type_col = ref_cell_type_col,
                                  cell_types = cell_types,
                                  dual_only = FALSE,
                                  n_cell_types = NULL)

    # Sort cell types for consistent reference category
    cell_types <- sort(cell_types)

    # Store reference cell type (first alphabetically)
    reference_cell_type <- cell_types[1]

    # Set dependent variables
    dep_vars <- paste0("PC", pc_subset)

    # Case 1: Query data only
    if(is.null(reference_data)){

        # Get query PCA variance for plotting
        query_pca_var <- attr(reducedDim(query_data, "PCA"), "percentVar")

        # Downsample reference and query data
        query_data <- downsampleSCE(sce_object = query_data,
                                    max_cells = max_cells_query,
                                    cell_types = cell_types,
                                    cell_type_col = query_cell_type_col)

        # Case 1a: No batch - PC ~ cell_type
        if(is.null(query_batch_col)){

            regress_data <- data.frame(
                reducedDim(query_data, "PCA")[, dep_vars],
                cell_type = factor(query_data[[query_cell_type_col]],
                                   levels = cell_types))

            # Regress PCs against cell types
            summaries <- lapply(dep_vars, regressFastCustom,
                                indep_var = "cell_type",
                                df = regress_data)
            names(summaries) <- dep_vars

            # Decompose R² by components
            r2_decomposed <- lapply(dep_vars, decomposeR2,
                                    indep_var = "cell_type",
                                    df = regress_data)
            names(r2_decomposed) <- dep_vars

            # Calculate variance contributions
            r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")
            var_expl <- query_pca_var[pc_subset]
            var_contr <- var_expl * r_squared
            total_var_expl <- sum(var_contr)

            # Calculate component-specific variance contributions
            r2_cell_type <- vapply(r2_decomposed, function(x) x$cell_type, numeric(1))
            var_contr_cell_type <- var_expl * r2_cell_type

            regress_res <- list(regression_summaries = summaries,
                                r_squared = r_squared,
                                r_squared_components = list(cell_type = r2_cell_type),
                                var_contributions = var_contr,
                                var_contributions_components = list(cell_type = var_contr_cell_type),
                                total_variance_explained = total_var_expl,
                                query_pca_var = query_pca_var,
                                indep_var = "cell_type",
                                reference_cell_type = reference_cell_type)

            regress_res <- adjustPValues(regress_res,
                                         adjust_method = adjust_method,
                                         indep_var = "cell_type")

            # Case 1b: With batch - PC ~ cell_type * batch
        } else {


            regress_data <- data.frame(
                reducedDim(query_data, "PCA")[, dep_vars],
                cell_type = factor(query_data[[query_cell_type_col]],
                                   levels = cell_types),
                batch = factor(query_data[[query_batch_col]]))

            # Get reference batch (first alphabetically)
            batch_levels <- levels(regress_data[["batch"]])
            reference_batch <- batch_levels[1]

            # Regress PCs against cell_type * batch interaction
            summaries <- lapply(dep_vars, regressFastCustom,
                                indep_var = "cell_type * batch",
                                df = regress_data)
            names(summaries) <- dep_vars

            # Decompose R² by components
            r2_decomposed <- lapply(dep_vars, decomposeR2,
                                    indep_var = "cell_type * batch",
                                    df = regress_data)
            names(r2_decomposed) <- dep_vars

            # Calculate variance contributions
            r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")
            var_expl <- query_pca_var[pc_subset]
            var_contr <- var_expl * r_squared
            total_var_expl <- sum(var_contr)

            # Calculate component-specific variance contributions
            r2_cell_type <- vapply(r2_decomposed, function(x) x$cell_type, numeric(1))
            r2_batch <- vapply(r2_decomposed, function(x) x$batch, numeric(1))
            r2_interaction <- vapply(r2_decomposed, function(x) x$interaction, numeric(1))

            var_contr_cell_type <- var_expl * r2_cell_type
            var_contr_batch <- var_expl * r2_batch
            var_contr_interaction <- var_expl * r2_interaction

            regress_res <- list(regression_summaries = summaries,
                                r_squared = r_squared,
                                r_squared_components = list(cell_type = r2_cell_type,
                                                            batch = r2_batch,
                                                            interaction = r2_interaction),
                                var_contributions = var_contr,
                                var_contributions_components = list(cell_type = var_contr_cell_type,
                                                                    batch = var_contr_batch,
                                                                    interaction = var_contr_interaction),
                                total_variance_explained = total_var_expl,
                                query_pca_var = query_pca_var,
                                indep_var = "cell_type_batch_interaction",
                                reference_cell_type = reference_cell_type,
                                reference_batch = reference_batch)

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
                                 cell_types = cell_types,
                                 pc_subset = pc_subset,
                                 assay_name = assay_name,
                                 max_cells_ref = max_cells_ref,
                                 max_cells_query = max_cells_query)

        # Case 2a: No batch - PC ~ cell_type * dataset
        if(is.null(query_batch_col)){

            pca_output[["dataset"]] <- factor(pca_output[["dataset"]],
                                              levels = c("Reference", "Query"))
            pca_output[["cell_type"]] <- factor(pca_output[["cell_type"]],
                                                levels = cell_types)

            # Unified interaction model: PC ~ cell_type * dataset
            summaries <- lapply(dep_vars, regressFastCustom,
                                indep_var = "cell_type * dataset",
                                df = pca_output)
            names(summaries) <- dep_vars

            # Decompose R² by components
            r2_decomposed <- lapply(dep_vars, decomposeR2,
                                    indep_var = "cell_type * dataset",
                                    df = pca_output)
            names(r2_decomposed) <- dep_vars

            # Calculate R-squared and variance contributions
            r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")
            var_expl <- reference_pca_var[pc_subset]
            var_contr <- var_expl * r_squared
            total_var_expl <- sum(var_contr)

            # Calculate component-specific variance contributions
            r2_cell_type <- vapply(r2_decomposed, function(x) x$cell_type, numeric(1))
            r2_dataset <- vapply(r2_decomposed, function(x) x$dataset, numeric(1))
            r2_interaction <- vapply(r2_decomposed, function(x) x$interaction, numeric(1))

            var_contr_cell_type <- var_expl * r2_cell_type
            var_contr_dataset <- var_expl * r2_dataset
            var_contr_interaction <- var_expl * r2_interaction

            regress_res <- list(regression_summaries = summaries,
                                r_squared = r_squared,
                                r_squared_components = list(cell_type = r2_cell_type,
                                                            dataset = r2_dataset,
                                                            interaction = r2_interaction),
                                var_contributions = var_contr,
                                var_contributions_components = list(cell_type = var_contr_cell_type,
                                                                    dataset = var_contr_dataset,
                                                                    interaction = var_contr_interaction),
                                total_variance_explained = total_var_expl,
                                reference_pca_var = reference_pca_var,
                                indep_var = "cell_type_dataset_interaction",
                                reference_cell_type = reference_cell_type)

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
                summaries <- lapply(dep_vars, regressFastCustom,
                                    indep_var = "cell_type * dataset",
                                    df = pca_output)
                names(summaries) <- dep_vars

                # Decompose R² by components
                r2_decomposed <- lapply(dep_vars, decomposeR2,
                                        indep_var = "cell_type * dataset",
                                        df = pca_output)
                names(r2_decomposed) <- dep_vars

                # Calculate R-squared and variance contributions
                r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")
                var_expl <- reference_pca_var[pc_subset]
                var_contr <- var_expl * r_squared
                total_var_expl <- sum(var_contr)

                # Calculate component-specific variance contributions
                r2_cell_type <- vapply(r2_decomposed, function(x) x$cell_type, numeric(1))
                r2_dataset <- vapply(r2_decomposed, function(x) x$dataset, numeric(1))
                r2_interaction <- vapply(r2_decomposed, function(x) x$interaction, numeric(1))

                var_contr_cell_type <- var_expl * r2_cell_type
                var_contr_dataset <- var_expl * r2_dataset
                var_contr_interaction <- var_expl * r2_interaction

                regress_res <- list(regression_summaries = summaries,
                                    r_squared = r_squared,
                                    r_squared_components = list(cell_type = r2_cell_type,
                                                                dataset = r2_dataset,
                                                                interaction = r2_interaction),
                                    var_contributions = var_contr,
                                    var_contributions_components = list(cell_type = var_contr_cell_type,
                                                                        dataset = var_contr_dataset,
                                                                        interaction = var_contr_interaction),
                                    total_variance_explained = total_var_expl,
                                    reference_pca_var = reference_pca_var,
                                    indep_var = "cell_type_dataset_interaction",
                                    reference_cell_type = reference_cell_type)

                regress_res <- adjustPValues(regress_res,
                                             adjust_method = adjust_method,
                                             indep_var = "cell_type_dataset_interaction")
            } else {
                # True multi-batch case: PC ~ cell_type * batch
                summaries <- lapply(dep_vars, regressFastCustom,
                                    indep_var = "cell_type * batch",
                                    df = pca_output)
                names(summaries) <- dep_vars

                # Decompose R² by components
                r2_decomposed <- lapply(dep_vars, decomposeR2,
                                        indep_var = "cell_type * batch",
                                        df = pca_output)
                names(r2_decomposed) <- dep_vars

                # Calculate R-squared and variance contributions
                r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")
                var_expl <- reference_pca_var[pc_subset]
                var_contr <- var_expl * r_squared
                total_var_expl <- sum(var_contr)

                # Calculate component-specific variance contributions
                r2_cell_type <- vapply(r2_decomposed, function(x) x$cell_type, numeric(1))
                r2_batch <- vapply(r2_decomposed, function(x) x$batch, numeric(1))
                r2_interaction <- vapply(r2_decomposed, function(x) x$interaction, numeric(1))

                var_contr_cell_type <- var_expl * r2_cell_type
                var_contr_batch <- var_expl * r2_batch
                var_contr_interaction <- var_expl * r2_interaction

                regress_res <- list(regression_summaries = summaries,
                                    r_squared = r_squared,
                                    r_squared_components = list(cell_type = r2_cell_type,
                                                                batch = r2_batch,
                                                                interaction = r2_interaction),
                                    var_contributions = var_contr,
                                    var_contributions_components = list(cell_type = var_contr_cell_type,
                                                                        batch = var_contr_batch,
                                                                        interaction = var_contr_interaction),
                                    total_variance_explained = total_var_expl,
                                    reference_pca_var = reference_pca_var,
                                    indep_var = "cell_type_batch_interaction",
                                    reference_cell_type = reference_cell_type)

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

#' @title Fast Custom Linear Regression for Principal Components
#'
#' @description
#' Performs efficient linear regression of principal component scores against categorical
#' predictors (cell types, batches, datasets, or their interactions) using QR decomposition
#' for numerical stability and computational efficiency.
#'
#' @details
#' This function implements a custom linear regression optimized for categorical predictors
#' commonly used in single-cell RNA sequencing analysis. It uses QR decomposition instead
#' of normal equations for improved numerical stability and handles rank-deficient design
#' matrices gracefully. The function supports various model specifications including:
#' \itemize{
#'   \item Simple cell type effects: \code{PC ~ cell_type}
#'   \item Cell type and batch interactions: \code{PC ~ cell_type * batch}
#'   \item Cell type and dataset interactions: \code{PC ~ cell_type * dataset}
#' }
#'
#' The output format is compatible with \code{speedglm} results to maintain consistency
#' with existing plotting and analysis workflows.
#'
#' @param pc A character string specifying the principal component column name in the data frame.
#' @param indep_var A character string specifying the independent variable specification.
#'   Options include "cell_type", "cell_type * batch", "cell_type * dataset", or other
#'   interaction specifications.
#' @param df A data frame containing the principal component scores and categorical predictors.
#'   Must include columns for the specified PC and predictor variables.
#'
#' @keywords internal
#'
#' @return A list containing:
#'   \item{coefficients}{A data frame with columns \code{coef}, \code{se}, \code{t}, and \code{p.value}
#'     containing regression coefficients, standard errors, t-statistics, and p-values.}
#'   \item{r_squared}{A numeric value representing the R-squared of the model.}
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @importFrom stats model.matrix pt
#'
# Perform linear regression for each principal component
regressFastCustom <- function(pc, indep_var, df) {

    # Extract the response variable
    y <- df[[pc]]
    n <- length(y)

    # Create design matrix based on independent variable specification
    if (indep_var == "cell_type") {
        X <- model.matrix(~ cell_type, data = df)
    } else if (indep_var == "cell_type * batch") {
        X <- model.matrix(~ cell_type * batch, data = df)
    } else if (indep_var == "cell_type * dataset") {
        X <- model.matrix(~ cell_type * dataset, data = df)
    } else {
        # Handle other interaction cases generically
        formula_str <- paste("~", gsub("_", " * ", gsub("_interaction$", "", indep_var)))
        X <- model.matrix(as.formula(formula_str), data = df)
    }

    # Use QR decomposition for numerical stability
    qr_decomp <- qr(X)

    # Handle rank deficiency (shouldn't happen with your data, but for safety)
    if (qr_decomp$rank < ncol(X)) {
        X <- X[, qr_decomp$pivot[1:qr_decomp$rank], drop = FALSE]
        qr_decomp <- qr(X)
    }

    # Solve for coefficients using QR decomposition
    coefficients <- qr.coef(qr_decomp, y)

    # Handle any NAs (shouldn't happen but for completeness)
    if (any(is.na(coefficients))) {
        # Fall back to normal equations
        XtX <- crossprod(X)
        Xty <- crossprod(X, y)
        coefficients <- solve(XtX, Xty)
    }

    # Compute fitted values and residuals
    fitted_values <- X %*% coefficients
    residuals <- y - fitted_values

    # Sum of squares calculations
    y_mean <- mean(y)
    ss_total <- sum((y - y_mean)^2)
    ss_residual <- sum(residuals^2)
    ss_regression <- ss_total - ss_residual

    # R-squared
    r_squared <- ss_regression / ss_total

    # Degrees of freedom
    p <- ncol(X)  # number of parameters including intercept
    df_residual <- n - p

    # Residual standard error
    mse <- ss_residual / df_residual

    # Variance-covariance matrix of coefficients
    # Var(beta) = sigma^2 * (X'X)^(-1)
    XtX_inv <- chol2inv(qr.R(qr_decomp))
    var_coef <- mse * XtX_inv
    se_coef <- sqrt(diag(var_coef))

    # t-statistics and p-values
    t_values <- coefficients / se_coef
    p_values <- 2 * pt(abs(t_values), df_residual, lower.tail = FALSE)

    # Create coefficients data.frame in exact format as speedlm
    coef_df <- data.frame(
        coef = as.numeric(coefficients),
        se = se_coef,
        t = t_values,
        p.value = p_values,
        stringsAsFactors = FALSE
    )

    # Set row names to match coefficient names
    rownames(coef_df) <- names(coefficients)

    # Return in same format as speedlm summary
    model_summary <- list(
        coefficients = coef_df,
        r_squared = r_squared
    )

    return(model_summary)
}

#' @title Decompose R-squared by Model Components
#'
#' @description
#' Decomposes the total R-squared from a linear model into individual components
#' representing the variance explained by main effects (cell type, batch/dataset)
#' and their interaction using sequential sum of squares.
#'
#' @details
#' This function performs R-squared decomposition using a sequential sum of squares approach,
#' which partitions the total explained variance into additive components. The decomposition
#' follows the hierarchical structure:
#' \enumerate{
#'   \item Cell type main effect
#'   \item Batch/dataset main effect (after accounting for cell type)
#'   \item Cell type × batch/dataset interaction (after accounting for main effects)
#' }
#'
#' For simple cell type models, only the cell type component is returned. For interaction
#' models, all three components are computed. The method uses group means and residual
#' analysis to avoid computationally expensive matrix operations while maintaining
#' mathematical accuracy equivalent to ANOVA decomposition.
#'
#' @param pc A character string specifying the principal component column name in the data frame.
#' @param indep_var A character string specifying the independent variable specification.
#'   Options are "cell_type", "cell_type * batch", or "cell_type * dataset".
#' @param df A data frame containing the principal component scores and categorical predictors.
#'   Must include columns for the specified PC and predictor variables.
#'
#' @keywords internal
#'
#' @return A named list containing R-squared components:
#'   \item{cell_type}{Numeric value representing the R-squared explained by cell type main effect.}
#'   \item{batch/dataset}{Numeric value representing the R-squared explained by batch or dataset
#'     main effect (only for interaction models).}
#'   \item{interaction}{Numeric value representing the R-squared explained by the interaction term
#'     (only for interaction models).}
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to decompose R-squared using sequential group means
decomposeR2 <- function(pc, indep_var, df) {
    y <- df[[pc]]
    y_mean <- mean(y)
    ss_total <- sum((y - y_mean)^2)

    if (indep_var == "cell_type") {
        # Simple case
        group_means <- tapply(y, df$cell_type, mean)
        y_pred <- group_means[df$cell_type]
        ss_explained <- sum((y_pred - y_mean)^2)

        return(list(cell_type = ss_explained / ss_total))

    } else if (indep_var %in% c("cell_type * batch", "cell_type * dataset")) {
        second_var <- if (indep_var == "cell_type * batch") "batch" else "dataset"

        # Sequential sum of squares approach
        # 1. Cell type only
        ct_means <- tapply(y, df$cell_type, mean)
        y_pred_ct <- ct_means[df$cell_type]
        ss_cell_type <- sum((y_pred_ct - y_mean)^2)

        # 2. Add second variable
        # Fit additive model by computing residuals and group means
        residuals_1 <- y - y_pred_ct
        sv_residual_means <- tapply(residuals_1, df[[second_var]], mean)
        y_pred_sv_adj <- sv_residual_means[df[[second_var]]]
        ss_second_var <- sum(y_pred_sv_adj^2)

        # 3. Add interaction (full model with all combinations)
        interaction_means <- tapply(y, list(df$cell_type, df[[second_var]]), mean)

        # Predict using full interaction model
        y_pred_full <- numeric(length(y))
        for (i in seq_along(y)) {
            ct_level <- as.character(df$cell_type[i])
            sv_level <- as.character(df[[second_var]][i])
            y_pred_full[i] <- interaction_means[ct_level, sv_level]
        }

        ss_full <- sum((y_pred_full - y_mean)^2)
        ss_interaction <- ss_full - ss_cell_type - ss_second_var

        result <- list(
            cell_type = ss_cell_type / ss_total,
            interaction = max(0, ss_interaction / ss_total)  # Ensure non-negative due to numerical precision
        )
        result[[second_var]] <- ss_second_var / ss_total

        return(result)
    }

    # Fallback
    return(list(total = 0))
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
