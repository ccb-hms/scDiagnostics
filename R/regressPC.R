#' Perform Linear Regression Analysis on Single-Cell Data
#'
#' This function performs linear regression on a SingleCellExperiment object, where the dependent variable is specified by the user as one of the principal components (PC1, PC2, etc.) from the dimension reduction slot, and the independent variable is provided as a column name in the SingleCellExperiment object.
#'
#' @param se_object A SingleCellExperiment object containing the data for regression analysis
#' @param dependent_vars A character string specifying the dependent variable principal component (e.g., "PC1", "PC2", etc.)
#' @param independent_var A character string specifying the column name for the independent variable in the SingleCellExperiment object
#'
#' @importFrom stats lm
#' @import SingleCellExperiment
#' @export
#'
#' @return The summary of the linear regression model
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(SingleR)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Run PCA
#' query_data <- runPCA(query_data)
#'
#' # Get cell type scoresx using SingleR or any other method
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#' 
#' # Specify the dependent variables (principal components) and independent variable (e.g., "labels")
#' dependent_vars <- c("PC1", "PC2", "PC3")
#' independent_var <- "labels"
#' 
#' # Perform linear regression on multiple principal components
#' result <- regressPC(se_object = query_data, 
#'                     dependent_vars = dependent_vars, 
#'                     independent_var = independent_var)
#' 
#' # Print the summaries of the linear regression models and R-squared values
#' print(result$regression_summaries)
#' print(result$rsquared_df)
#' 
#' # Note: Instead of using SingleR, you can use any other method 
#' # to obtain the scores for regression analysis.
#'
regressPC <- function(se_object, 
                      dependent_vars, 
                      independent_var) {
  
  # Get the dependent variables from the specified principal components
  dependent_list <- lapply(dependent_vars, function(var) {
    reducedDim(se_object, "PCA")[, var, drop = FALSE]
  })
  
  # Create a data frame with the dependent and independent variables
  df <- data.frame(Independent = colData(se_object)[[independent_var]])
  for (i in seq_along(dependent_list)) {
    df[[paste0("PC", i)]] <- dependent_list[[i]]
  }
  
  # Perform linear regression for each principal component
  regression_summaries <- list()
  for (i in seq_along(dependent_list)) {
    lm_model <- lm(paste0("PC", i, " ~ Independent"), data = df)
    regression_summary <- summary(lm_model)
    regression_summaries[[paste0("PC", i)]] <- regression_summary
  }
  
  # Calculate R-squared values
  rsquared <- sapply(regression_summaries, function(summary) {
    rsq <- summary$r.squared
    return(rsq)
  })
  
  rsquared_df <- data.frame(PC = names(regression_summaries), R2 = rsquared)
  
  # Return both the summaries of the linear regression models and R-squared values
  return(list(regression_summaries = regression_summaries, rsquared_df = rsquared_df))
}

