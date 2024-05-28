#' @title Plot Visualization of Output from comparePCASubspace Function
#' 
#' @description This function generates a visualization of the output from the `comparePCASubspace` function.
#' The plot shows the cosine of principal angles between reference and query principal components,
#' with point sizes representing the variance explained.
#' 
#' @details The function converts the input list into a data frame suitable for plotting with `ggplot2`.
#' Each point in the scatter plot represents the cosine of a principal angle, with the size of the point
#' indicating the average variance explained by the corresponding principal components.
#' 
#' @param x A numeric matrix output from the `comparePCA` function, representing 
#' cosine similarities between query and reference principal components.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot object representing the heatmap of cosine similarities.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{comparePCASubspace}}
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
#' # Compare PCA subspaces
#' subspace_comparison <- comparePCASubspace(query_data_subset, ref_data_subset, 
#'                                           pc_subset = c(1:5))
#' 
#' # Create a data frame for plotting
#' plot(subspace_comparison)
#' 
# Function to produce the visualization of output from comparePCASubspace function
plot.comparePCASubspace <- function(x, ...){
    
    # Create a data frame for plotting
    x <- data.frame(PC = paste0("Ref PC", subspace_comparison$cosine_id[, 1],
                                              " - Query PC", subspace_comparison$cosine_id[, 2]),
                                  Cosine = subspace_comparison$cosine_similarity,
                                  VarianceExplained = subspace_comparison$var_explained_avg)
    x$PC <- factor(x$PC, levels = x$PC)
    
    # Create plot
    pc_plot <- ggplot2::ggplot(x, aes(x = PC, y = Cosine, size = VarianceExplained)) +
        ggplot2::geom_point() +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(title = "Principal Angles Cosines with Variance Explained",
                      x = "",
                      y = "Cosine of Principal Angle",
                      size = "Variance Explained") +
        ggplot2::theme_minimal()
    print(pc_plot)
}