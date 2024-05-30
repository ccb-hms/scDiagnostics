#' @title 
#' Plot the output of the calculateAveragePairwiseCorrelation function
#'
#' @description 
#' This function takes the output of the calculateAveragePairwiseCorrelation function,
#' which should be a matrix of pairwise correlations, and plots it as a heatmap.
#' 
#' @details 
#' This function converts the correlation matrix into a dataframe, creates a heatmap using ggplot2,
#' and customizes the appearance of the heatmap with updated colors and improved aesthetics.
#'
#' @param x Output matrix from calculateAveragePairwiseCorrelation function.
#' @param ... Additional arguments to be passed to the plotting function.
#'
#' @return A ggplot2 object representing the heatmap plot.
#' 
#' @export
#'         
#' @seealso \code{\link{calculateAveragePairwiseCorrelation}}
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
#' set.seed(100)
#' indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
#' ref_data <- sce[, indices]
#' query_data <- sce[, -indices]
#'
#' # log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR
#' scores <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- scores$labels
#'
#' # Compute Pairwise Correlations
#' # Note: The selection of highly variable genes and desired cell types may vary 
#' # based on user preference. 
#' # The cell type annotation method used in this example is SingleR. 
#' # User can use any other method for cell type annotation and provide 
#' # the corresponding labels in the metadata.
#'
#' # Selecting highly variable genes
#' ref_var <- getTopHVGs(ref_data, n = 2000)
#' query_var <- getTopHVGs(query_data, n = 2000)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#'
#' # Select desired cell types
#' selected_cell_types <- c("CD4", "CD8", "B_and_plasma")
#' ref_data_subset <- ref_data[common_genes, ref_data$reclustered.broad %in% selected_cell_types]
#' query_data_subset <- query_data[common_genes, query_data$reclustered.broad %in% selected_cell_types]
#' 
#' # Run PCA on the reference data
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' # Compute pairwise correlations
#' cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data_subset, 
#'                                                       reference_data = ref_data_subset, 
#'                                                       n_components = 10,
#'                                                       query_cell_type_col = "labels", 
#'                                                       ref_cell_type_col = "reclustered.broad", 
#'                                                       cell_types = selected_cell_types, 
#'                                                       correlation_method = "spearman")
#'
#' # Visualize the results
#' plot(cor_matrix_avg)
#' 
#'
# Function to plot the output of the calculateAveragePairwiseCorrelation function
plot.calculateAveragePairwiseCorrelation <- function(x, ...){
    
    # Convert matrix to dataframe
    cor_df <- as.data.frame(as.table(cor_matrix_avg))
    cor_df$Var1 <- factor(cor_df$Var1, levels = rownames(cor_matrix_avg))
    cor_df$Var2 <- factor(cor_df$Var2, levels = rev(colnames(cor_matrix_avg)))
    
    # Create the heatmap with updated colors and improved aesthetics
    heatmap_plot <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var2, y = Var1)) +
        ggplot2::geom_tile(ggplot2::aes(fill = Freq), color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = round(Freq, 2)), color = "black", size = 3, family = "sans") +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                                      midpoint = 0, limits = c(min(cor_df$Freq), max(cor_df$Freq)),
                                      name = "Correlation",
                                      breaks = seq(-1, 1, by = 0.2)) +  # Specify color scale breaks
        ggplot2::labs(title = "Correlation Heatmap", x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
                       axis.text.y = ggplot2::element_text(family = "sans"),  # Set font family for y-axis labels
                       plot.title = ggplot2::element_text(face = "bold"),  # Make title bold
                       legend.position = "right",  # Place legend on RHS
                       legend.title = ggplot2::element_text(face = "italic"))
    
    # Print the plot
    print(heatmap_plot)
}
