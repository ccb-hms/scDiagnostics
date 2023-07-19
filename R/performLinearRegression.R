#' Perform Linear Regression Analysis on Single-Cell Data
#'
#' This function performs linear regression on a reaference or query object, where the dependent variable is specified by the user as one of the principal components (PC1, PC2, etc.) from the dimension reduction slot, and the independent variable is provided as a column name in the colData of SingleCellExperiment object.
#'
#' @param se_object A SingleCellExperiment object containing the data for regression analysis
#' @param dependent_var A character string specifying the dependent variable principal component (e.g., "PC1", "PC2", etc.)
#' @param independent_var A character string specifying the column name for the independent variable in the colData of SingleCellExperiment object
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
#' colData(query_data)$labels <- pred$labels
#'
#' # Perform linear regression on PC1 and a specific variable (e.g., "labels")
#' performLinearRegression(query_data, "PC1", "labels")
#'
#' # Note: Instead of using SingleR, you can use any other method to obtain the cell labels for regression analysis.
#'
performLinearRegression <- function(se_object, dependent_var, independent_var) {

  # Get the dependent variable from the specified principal component
  dependent <- reducedDim(se_object, "PCA")[, dependent_var]

  # Create a data frame with the dependent and independent variables
  df <- data.frame(Dependent = dependent,
                   Independent = colData(se_object)[[independent_var]])

  # Perform linear regression
  lm_model <- lm(Dependent ~ Independent, data = df)

  # Print and return the summary of the linear regression model
  regression_summary <- summary(lm_model)
  print(regression_summary)

  return(regression_summary)
}
