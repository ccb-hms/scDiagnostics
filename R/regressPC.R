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
#' Briefly, the \eqn{R^2} is calculated from a linear regression of the covariate \eqn{B} of 
#' interest onto each principal component. The variance contribution of the 
#' covariate effect per principal component is then calculated as the product of 
#' the variance explained by the i-th principal component (PC) and the 
#' corresponding \eqn{R^2(PC_i | B)}. The sum across all variance contributions by the 
#' covariate effects in all principal components gives the total variance 
#' explained by the covariate as follows:
#' 
#' \deqn{Var(C|B) = \sum_{i=1}^G \text{Var}(C|PC_i) \times R^2(PC_i | B)}
#' 
#' where, \eqn{\text{Var}(C \mid PC_i)} is the variance of the data matrix \eqn{C} 
#' explained by the i-th principal component. See references for details.
#' 
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' If NULL, the PC scores are regressed against the cell types of the reference data.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#' @param adjust_method A character string specifying the method to adjust the p-values. 
#'   Options include "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", or "none". 
#'   Default is "BH" (Benjamini-Hochberg). Default is "BH".
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
#' Andrew Ghazi, \email{andrew_ghazi@hms.harvard.edu}
#' 
#' @seealso \code{\link{plot.regressPC}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Plot the PC data (no query data)
#' regress_res <- regressPC(reference_data = reference_data,
#'                          ref_cell_type_col = "expert_annotation", 
#'                          cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                          pc_subset = 1:15)
#' # Plot results
#' plot(regress_res, plot_type = "r_squared")
#' plot(regress_res, plot_type = "p-value")
#'
#' # Plot the PC data (with query data)
#' regress_res <- regressPC(reference_data = reference_data,
#'                          query_data = query_data,
#'                          ref_cell_type_col = "expert_annotation", 
#'                          query_cell_type_col = "SingleR_annotation",
#'                          cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                          pc_subset = 1:15)
#' # Plot results
#' plot(regress_res, plot_type = "r_squared")
#' plot(regress_res, plot_type = "p-value")
#'
#' @importFrom rlang .data
#' @import SingleCellExperiment
#' 
regressPC <- function(reference_data,
                      query_data = NULL, 
                      ref_cell_type_col, 
                      query_cell_type_col = NULL, 
                      cell_types = NULL,
                      pc_subset = 1:10,
                      adjust_method = c("BH", "holm", 
                                        "hochberg", "hommel", 
                                        "bonferroni", "BY", 
                                        "fdr", "none")) {
    
    # Match argument for independent variable
    adjust_method <- match.arg(adjust_method)
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset)
    
    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        if(is.null(query_data)){
            cell_types <- na.omit(unique(
                c(reference_data[[ref_cell_type_col]])))
        } else{
            cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                           query_data[[query_cell_type_col]])))
        }
    }
    
    # Perform linear regression for each principal component
    .regressFast <- function(pc, indep_var, df, reference_category) {
        regression_formula <- paste(pc, "~", indep_var)
        if(isTRUE(reference_category)){
            model <- do.call(speedglm::speedlm, 
                             list(formula = paste(pc, " ~ ", indep_var), 
                                  data = df))
        } else{
            model <- do.call(speedglm::speedlm, 
                             list(formula = paste(pc, " ~ ", indep_var, " - 1"), 
                                  data = df))
        }
        model_summary <- list(coefficients = summary(model)[["coefficients"]],
                              r_squared = summary(model)[["r.squared"]])
        return(model_summary)
    }
    
    # Set dependent variables
    dep_vars <- paste0("PC", pc_subset)
    
    # Regress PCs against cell types
    if(is.null(query_data)){
        
        # Create a data frame with the dependent and independent variables
        reference_labels <- reference_data[[ref_cell_type_col]]
        regress_data <- data.frame(
            reducedDim(reference_data, "PCA")[, dep_vars], 
            Cell_Type_ = reference_data[[ref_cell_type_col]])
        regress_data <- regress_data[reference_labels %in% cell_types,]
        
        # Regress PCs
        summaries <- lapply(dep_vars, .regressFast, indep_var = "Cell_Type_", 
                            df = regress_data, 
                            reference_category = TRUE)
        names(summaries) <- dep_vars
        
        # Adjust cell types for the summaries
        for(pc in seq_len(length(summaries))){
            rownames(summaries[[pc]][["coefficients"]]) <- sort(cell_types)
        }
        
        # Calculate variance contributions by principal component
        r_squared <- vapply(summaries, `[[`, numeric(1), x = "r_squared")
        var_expl <- attr(reducedDim(reference_data, "PCA"), 
                         "percentVar")[pc_subset]
        var_contr <- var_expl * r_squared
        
        # Calculate total variance explained by summing the variance contributions
        total_var_expl <- sum(var_contr)
        
        # Return the summaries of the linear regression models, R-squared values, and variance contributions
        regress_res <- list(regression_summaries = summaries,
                            r_squared = r_squared,
                            var_contributions = var_contr,
                            total_variance_explained = total_var_expl,
                            indep_var = "cell_type")
        
        # Adjust p-values
        regress_res <- adjustPValues(regress_res, 
                                     adjust_method = adjust_method, 
                                     indep_var = "cell_type")
        
    } else {
        
        # Get the projected PCA data
        pca_output <- projectPCA(query_data = query_data, 
                                 reference_data = reference_data, 
                                 query_cell_type_col = query_cell_type_col, 
                                 ref_cell_type_col = ref_cell_type_col,
                                 pc_subset = pc_subset)
        pca_output <- pca_output[pca_output[["cell_type"]] %in% cell_types,]
        pca_output[["dataset"]] <- factor(pca_output[["dataset"]], 
                                          levels = c("Reference", "Query"))
        
        # Set dependent variables
        dep_vars <- paste0("PC", pc_subset)
        
        # Set independent variables
        indep_list <- split(pca_output, pca_output[["cell_type"]])
        
        # Adding regression summaries for each cell type
        regress_res <- vector("list", length = length(indep_list))
        names(regress_res) <- cell_types
        for(cell_type in cell_types){
            
            regress_res[[cell_type]] <- lapply(dep_vars, .regressFast, 
                                               indep_var = "dataset", 
                                               df = indep_list[[cell_type]],
                                               reference_category = TRUE)
            names(regress_res[[cell_type]]) <- dep_vars
        }
        regress_res[["indep_var"]] <- "dataset"
        
        # Adjust p-values
        regress_res <- adjustPValues(regress_res, 
                                     adjust_method = adjust_method,
                                     indep_var = "dataset")
    }
    
    # Return regression output
    class(regress_res) <- c(class(regress_res), "regressPC")
    return(regress_res)
}

#' @title Adjust P-Values in Regression Results
#'
#' @description
#' Adjusts the p-values in the regression results using a specified adjustment method.
#' The adjustment is performed either for each principal component (PC) by cell type 
#' or for each dataset, depending on the selected independent variable.
#'
#' @details
#' This function adjusts p-values from regression results stored in a list. The adjustment 
#' can be applied either across cell types or datasets, depending on the user’s choice. 
#' The method for adjusting p-values can be selected from various options such as Benjamini-Hochberg (BH),
#' Holm, and others, which are supported by the `p.adjust` function in R.
#'
#' @param regress_res A list containing regression results. The structure of the list 
#'   depends on the \code{indep_var} argument: if \code{indep_var} is "cell_type", 
#'   the list should contain regression summaries for each principal component (PC);
#'   if \code{indep_var} is "dataset", it should contain summaries for each dataset.
#' @param adjust_method A character string specifying the method to adjust the p-values. 
#'   Options include "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", or "none". 
#'   Default is "BH" (Benjamini-Hochberg). Default is "BH".
#' @param indep_var A character string specifying the independent variable for the adjustment. 
#'   Options are "cell_type" (default) or "dataset".
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
# Function to adjust p-values from the regression analysis on PCs
adjustPValues <- function(regress_res, 
                          adjust_method = c("BH", "holm", 
                                            "hochberg", "hommel", 
                                            "bonferroni", "BY", 
                                            "fdr", "none"),
                          indep_var = c("cell_type", "dataset")){
    
    # Match argument for independent variable
    adjust_method <- match.arg(adjust_method)
    
    # Match argument for independent variable
    indep_var <- match.arg(indep_var)
    
    if(indep_var == "cell_type"){
        
        # Add adjusted p-values 
        for (pc in names(regress_res$regression_summaries)) {
            
            # Extract the coefficients table for the current PC
            coeffs <- regress_res$regression_summaries[[pc]]$coefficients
            
            # Adjust the p-values using Benjamini-Hochberg (FDR) correction
            p.adjusted <- p.adjust(coeffs$p.value, method = "BH")
            
            # Add the adjusted p-values as a new column to the coefficients table
            coeffs$p.adjusted <- p.adjusted
            
            # Update the coefficients table in the object
            regress_res$regression_summaries[[pc]]$coefficients <- coeffs
        }
        
    } else if (indep_var == "dataset"){
        
        # Add adjusted p-values 
        for (cell_type in names(regress_res)[-length(regress_res)]) {
            
            for(pc_id in seq_len(length(regress_res[[cell_type]]))){
                
                # Extract the coefficients table for the current PC
                coeffs <- regress_res[[cell_type]][[pc_id]]$coefficients
                
                # Adjust the p-values using Benjamini-Hochberg (FDR) correction
                p.adjusted <- p.adjust(coeffs$p.value, method = "BH")
                
                # Add the adjusted p-values as a new column to the coefficients table
                coeffs$p.adjusted <- p.adjusted
                
                # Update the coefficients table in the object
                regress_res[[cell_type]][[pc_id]]$coefficients <- coeffs
            }
        }
    }
    
    # Return regression object
    return(regress_res)
}
