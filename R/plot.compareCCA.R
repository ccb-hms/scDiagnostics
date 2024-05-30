#' @title Plot Visualization of Output from compareCCA Function
#' 
#' @description This function generates a visualization of the output from the `compareCCA` function.
#' The plot shows the cosine similarities of canonical correlation analysis (CCA) coefficients,
#' with point sizes representing the correlations.
#'
#' @details The function converts the input list into a data frame suitable for plotting with `ggplot2`.
#' Each point in the scatter plot represents the cosine similarity of CCA coefficients, with the size of the point
#' indicating the correlation.
#'
#' @param x A list containing the output from the `compareCCA` function. 
#' This list should include `cosine_similarity` and `correlations`.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the scatter plot of cosine similarities of CCA coefficients and correlations.
#'
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{compareCCA}}
#' 
#' @examples
#' # Load necessary library
#' library(scRNAseq)
#' library(scuttle)
#' library(scran)
#' library(SingleR)
#' library(ggplot2)
#' library(scater)
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
#' # Compare CCA
#' cca_comparison <- compareCCA(query_data_subset, ref_data_subset, 
#'                              pc_subset = c(1:5))
#' 
#' # Visualize output of CCA comparison
#' plot(cca_comparison)
#' 
#' 
# Plot visualization of output from compareCCA function
plot.compareCCA <- function(x, ...){
    
    # Create a data frame for plotting
    comparison_data <- data.frame(CCA = paste0("CC", 1:length(x$correlations)),
                                  Cosine = x$cosine_similarity,
                                  Correlation = x$correlations)
    comparison_data$CC <- factor(comparison_data$CCA, levels = comparison_data$CCA)
    
    
    cca_plot <- ggplot2::ggplot(comparison_data, aes(x = CCA, y = Cosine, size = Correlation)) +
        ggplot2::geom_point() +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(title = "Cosine Similarities of CCA Coefficients with Correlation",
                      x = "",
                      y = "Cosine of CC Coefficients",
                      size = "Correlation") +
        ggplot2::theme_minimal()
    print(cca_plot)
}