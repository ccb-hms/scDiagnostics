#' Plot gene expression distribution from overall and cell type-specific perspective
#'
#' This function generates histogram plots to visualize the distribution of gene expression values for a specific gene
#' across the overall dataset and within a specified cell type.
#'
#' @param reference_data An object of class "SingleCellExperiment" containing log-transformed expression matrix and metadata for the reference dataset.
#' @param query_data An object of class "SingleCellExperiment" containing log-transformed expression matrix and metadata for the query dataset.
#' @param reference_cell_labels A character column name from the colData(reference_data) of the reference dataset representing cell type information.
#' @param query_cell_labels A character column name from the colData(query_data) of the query dataset representing cell type information.
#' @param gene_name A character string representing the gene name for which the distribution is to be visualized.
#' @param label A character vector of cell type labels to plot (e.g., c("T-cell", "B-cell"))

#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assay
#' @import SingleCellExperiment
#'
#' @return A object displaying density plots of gene expression distribution from both the 
#'         overall dataset and a cell type-specific perspective.
#'         
#' @export
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
#' # Log transform datasets
#' ref_data <- logNormCounts(ref_data)
#' query_data <- logNormCounts(query_data)
#'
#' # Get cell type scores using SingleR or any other method
#' pred <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
#'
#' # Add labels to query object
#' colData(query_data)$labels <- pred$labels
#'
#' # Note: Users can use SingleR or any other method to obtain the cell type annotations.
#' plotMarkerExpression(reference_data = ref_data, 
#'                      query_data = query_data, 
#'                      reference_cell_labels = "reclustered.broad", 
#'                      query_cell_labels = "labels", 
#'                      gene_name = "VPREB3", 
#'                      label = "B_and_plasma")
#'
plotMarkerExpression <- function(reference_data, 
                                 query_data, 
                                 reference_cell_labels, 
                                 query_cell_labels, 
                                 gene_name, 
                                 label) {
  
  # Get expression of the specified gene for reference and query datasets
  reference_gene_expression <- assay(reference_data, "logcounts")[gene_name, ]
  query_gene_expression <- assay(query_data, "logcounts")[gene_name, ]
  
  # Create a combined vector of gene expression values
  combined_gene_expression <- c(reference_gene_expression, query_gene_expression)
  
  # Create a grouping vector for dataset labels
  dataset_labels <- rep(c("Reference", "Query"), times = c(length(reference_gene_expression), length(query_gene_expression)))
  
  # Combine the gene expression values and dataset labels
  data <- data.frame(
    GeneExpression = combined_gene_expression,
    Dataset = dataset_labels
  )
  
  # Create a stacked density plot using ggplot2 for overall dataset
  overall_plot <- ggplot(data, aes(x = GeneExpression, fill = Dataset)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Overall Distribution"), x = paste("Log gene Expression", gene_name), y = "Density") +
    theme_minimal()
  
  # Create a subset of data for cell type-specific distribution
  index1 <- which(reference_data[[reference_cell_labels]] %in%  label)
  index2 <- which(query_data[[query_cell_labels]] %in%  label)
  
  reference_gene_expression_cell_type <- assay(reference_data, "logcounts")[gene_name, index1]
  query_gene_expression_cell_type <- assay(query_data, "logcounts")[gene_name, index2]
  
  # Combine the gene expression values and dataset labels for cell type-specific
  combined_gene_expression <- c(reference_gene_expression_cell_type, query_gene_expression_cell_type)
  
  # Create a grouping vector for dataset labels
  dataset_labels <- rep(c("Reference", "Query"), times = c(length(reference_gene_expression_cell_type), length(query_gene_expression_cell_type)))
  
  # Combine the gene expression values and dataset labels
  cell_type_specific_data <- data.frame(
    GeneExpression = combined_gene_expression,
    Dataset = dataset_labels
  )
  
  # Create a stacked density plot using ggplot2 for cell type-specific dataset
  cell_type_specific_plot <- ggplot(cell_type_specific_data, aes(x = GeneExpression, fill = Dataset)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Cell Type-Specific Distribution"), x = paste("Log gene Expression", gene_name), y = "Density") +
    theme_minimal()
  
  return(grid.arrange(overall_plot, cell_type_specific_plot, ncol = 2))
}