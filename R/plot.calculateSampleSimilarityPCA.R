#' @title Plot Cosine Similarities Between Samples and PCs
#'
#' @description 
#' This function creates a heatmap plot to visualize the cosine similarities between samples and principal components (PCs).
#'
#' @details 
#' This function reshapes the input data frame to create a long format suitable for plotting as a heatmap. It then
#' creates a heatmap plot using ggplot2, where the x-axis represents the PCs, the y-axis represents the samples, and the
#' color intensity represents the cosine similarity values.
#'
#' @param x An object of class 'calculateSampleSimilarityPCA' containing a dataframe of cosine similarity values 
#' between samples and PCs.
#' @param pc_subset A numeric vector specifying the subset of principal components to include in the plot (default: c(1:5)).
#'
#' @return A ggplot object representing the cosine similarity heatmap.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{calculateSampleSimilarityPCA}}
#'
#' @examples
#' # Load required libraries
#' library(scRNAseq)
#' library(scuttle)
#' library(SingleR)
#' library(scran)
#' library(scater)
#'
#' # Load data (replace with your data loading)
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
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
#' ref_var <- scran::getTopHVGs(ref_data, n = 2000)
#' query_var <- scran::getTopHVGs(query_data, n = 2000)
#' 
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#' ref_data_subset <- ref_data[common_genes, ]
#' query_data_subset <- query_data[common_genes, ]
#'
#' # Run PCA on the reference data (assumed to be prepared)
#' ref_data_subset <- runPCA(ref_data_subset)
#'
#' # Store PCA anomaly data and plots
#' anomaly_output <- detectAnomaly(reference_data = ref_data_subset, 
#'                                 ref_cell_type_col = "reclustered.broad", 
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5) 
#' top6_anomalies <- names(sort(anomaly_output$Combined$reference_anomaly_scores, decreasing = TRUE)[1:6])
#' 
#' # Compute cosine similarity between anomalies and top PCs
#' cosine_similarities <- calculateSampleSimilarityPCA(ref_data_subset, samples = top6_anomalies, 
#'                                                     pc_subset = c(1:10), n_top_vars = 50)
#' cosine_similarities
#' 
#' # Plot similarities
#' plot(cosine_similarities, pc_subset = c(1:5))
#' 
# Function to plot cosine similarities between samples and PCs
plot.calculateSampleSimilarityPCA <- function(x, pc_subset = c(1:5), ...){
    
    # Subset data
    x <- x[, paste0("PC", pc_subset)]
    
    # Initialize empty vectors for reshaped data
    sample_names <- c()
    pc_names <- c()
    cosine_values <- c()
    
    # Loop through the data frame to manually reshape it
    for (sample in rownames(x)) {
        for (pc in colnames(x)) {
            sample_names <- c(sample_names, sample)
            pc_names <- c(pc_names, pc)
            cosine_values <- c(cosine_values, x[sample, pc])
        }
    }
    
    # Create a data frame with the reshaped data
    cosine_long <- data.frame(Sample = factor(sample_names, levels = rev(rownames(x))), 
                              PC = pc_names, CosineSimilarity = cosine_values)
    
    # Create the heatmap plot
    plot <- ggplot(cosine_long, aes(x = PC, y = Sample, fill = CosineSimilarity)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.2f", CosineSimilarity)), size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                             limits = c(-1, 1), space = "Lab", name = "Cosine Similarity") +
        labs(title = "Cosine Similarity Heatmap",
             x = "",
             y = "") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5))
    return(plot)
}

