#' @title Plot Distance Density Comparison for a Specific Cell Type and Sample
#'
#' @description This function plots the density functions for the reference data and the distances from a specified query sample 
#' to all reference samples within a specified cell type.
#'
#' @details The function first checks if the specified cell type and sample name are present in the \code{distance_data}. If the 
#' specified cell type or sample name is not found, an error is thrown. It then extracts the distances within the reference dataset 
#' and the distances from the specified query sample to the reference samples. The function creates a density plot using \code{ggplot2} 
#' to compare the distance distributions. The density plot will show two distributions: one for the pairwise distances within the 
#' reference dataset and one for the distances from the specified query sample to each reference sample. These distributions are 
#' plotted in different colors to visually assess how similar the query sample is to the reference samples of the specified cell type.
#'
#' @param distance_data A list containing the distance data computed by \code{computeDistanceDiagnostics}.
#' @param cell_type A string specifying the cell type to plot.
#' @param sample_name A string specifying the query sample name for which to plot the distances.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A ggplot2 density plot comparing the reference distances and the distances from the specified sample to the reference samples.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
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
#' # Plot the PC data
#' distance_data <- computeDistanceDiagnostics(query_data, reference_data, 
#'                                             n_components = 10, 
#'                                             query_cell_type_col = "labels", 
#'                                             ref_cell_type_col = "reclustered.broad",
#'                                             pc_subset = c(1:10)) 
#' 
#' # Identify outliers for CD4
#' cd4_anomalites <- detectAnomaly(reference_data = ref_data_subset, query_data = query_data_subset, 
#'                                 query_cell_type_col = "labels", 
#'                                 ref_cell_type_col = "reclustered.broad",
#'                                 n_components = 10,
#'                                 n_tree = 500,
#'                                 anomaly_treshold = 0.5,
#'                                 store_plots = TRUE)$CD4
#' cd4_top_anomaly <- names(which.max(cd4_anomalites$anomaly_scores))
#' 
#' # Plot the densities of the distances
#' plot(distance_data, cell_type = "CD4", sample_name = cd4_top_anomaly)
#'  
# Function to plot density functions for the reference data and the specified sample
plot.computeDistanceDiagnostics <- function(distance_data, cell_type, sample_name, ...) {
    
    # Check if cell type is available
    if(length(cell_type) != 1 || !(cell_type %in% names(distance_data)))
        stop("The specified \'cell_type\' is not available.")
    
    # Filter distance data for the specified cell type
    cell_distances <- distance_data[[cell_type]]
    
    # Check if sample is available in data for that cell type
    if(!(sample_name %in% rownames(cell_distances$query_to_ref_distances)))
        stop("The specified \'cell_name\' is not available for that cell type.")
    
    # Extract distances within the reference dataset
    ref_distances <- cell_distances$ref_distances
    
    # Extract distances from the specified sample to reference samples
    sample_distances <- cell_distances$query_to_ref_distances[sample_name,]
    
    # Create data frame for plotting
    plot_data <- data.frame(Distance_Type = c(rep("Reference", length(ref_distances)), rep("Sample", length(sample_distances))),
                            Distance = c(ref_distances, sample_distances))
    
    # Plot density comparison
    density_plot <- ggplot(plot_data, aes(x = Distance, fill = Distance_Type)) +
        geom_density(alpha = 0.5) +
        labs(title = paste("Density Comparison for Cell Type:", cell_type, "and Sample:", sample_name),
             x = "Distance", y = "Density") +
        scale_fill_manual(name = "Distance Type", values = c("Reference" = "blue", "Sample" = "red")) +
        theme_minimal()
    
    # Print the density plot
    print(density_plot)
}







