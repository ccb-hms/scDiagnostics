#' @title Plot gene expression distribution from overall and cell type-specific perspective
#' 
#' @description
#' This function generates density plots to visualize the distribution of gene expression values 
#' for a specific gene across the overall dataset and within a specified cell type.
#'
#' @details 
#' This function generates density plots to compare the distribution of a specific marker 
#' gene between reference and query datasets. The aim is to inspect the alignment of gene expression 
#' levels as a surrogate for dataset similarity. Similar distributions suggest a good alignment, 
#' while differences may indicate discrepancies or incompatibilities between the datasets.
#' 
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param gene_name character. A string representing the gene name for which the distribution is to be visualized.
#' @param label character. A vector of cell type labels to plot (e.g., c("T-cell", "B-cell")).
#'
#' @return A gtable object containing two arranged density plots as grobs. 
#'         The first plot shows the overall gene expression distribution, 
#'         and the second plot displays the cell type-specific expression 
#'         distribution.
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
#'                      ref_cell_type_col = "reclustered.broad", 
#'                      query_cell_type_col = "labels", 
#'                      gene_name = "VPREB3", 
#'                      label = "B_and_plasma")
#' 
#'
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assay
#' @import SingleCellExperiment
#' @export
plotMarkerExpression <- function(reference_data, 
                                 query_data, 
                                 ref_cell_type_col, 
                                 query_cell_type_col, 
                                 gene_name, 
                                 label) {
  # Sanity checks
  # Check if query_data is a SingleCellExperiment object
  if (!is(query_data, "SingleCellExperiment")) {
    stop("query_data must be a SingleCellExperiment object.")
  }
  
  # Check if reference_data is a SingleCellExperiment object
  if (!is(reference_data, "SingleCellExperiment")) {
    stop("reference_data must be a SingleCellExperiment object.")
  }
  
  # Check if gene_name is present in both query_data and reference_data
  if (!(gene_name %in% rownames(assay(query_data)) && gene_name %in% 
        rownames(assay(reference_data)))) {
    stop("gene_name: '", gene_name, "' is not present in the 
         row names of both query_data and reference_data.")
  }
    
  # Check if all labels are present in query_data
  if (!all(label %in% query_data[[query_cell_type_col]])) {
    stop("One or more labels specified are not present in query_data.")
  }
  
  # Check if all labels are present in reference_data
  if (!all(label %in% reference_data[[ref_cell_type_col]])) {
    stop("One or more labels specified are not present in reference_data.")
  }
  
  # Get expression of the specified gene for reference and query datasets
  ref_gene_expression <- assay(reference_data, "logcounts")[gene_name, ]
  query_gene_expression <- assay(query_data, "logcounts")[gene_name, ]
  ref_gene_expression_specific <- assay(reference_data, "logcounts")[gene_name, which(reference_data[[ref_cell_type_col]] %in% label)]
  query_gene_expression_specific <- assay(query_data, "logcounts")[gene_name, which(query_data[[query_cell_type_col]] %in% label)]
  
  # Create a combined vector of gene expression values
  combined_gene_expression <- c(ref_gene_expression, query_gene_expression,
                                ref_gene_expression_specific, query_gene_expression_specific)
  
  # Create a grouping vector for dataset labels
  dataset_labels <- rep(c("Reference", "Query", "Reference", "Query"), 
                        times = c(length(ref_gene_expression), length(query_gene_expression),
                                  length(ref_gene_expression_specific), length(query_gene_expression_specific)))
  
  # Combine the gene expression values and dataset labels
  data <- data.frame(GeneExpression = combined_gene_expression,
                     Dataset = dataset_labels,
                     plot_type = rep(c("Overall Distribution", "Cell Type-Specific Distribution"), 
                                     times = c(length(ref_gene_expression) + length(query_gene_expression),
                                               length(ref_gene_expression_specific) + length(query_gene_expression_specific))))
  
  # Create a stacked density plot using ggplot2 for overall dataset
  plot_obj <- ggplot2::ggplot(data, ggplot2::aes(x = GeneExpression, fill = Dataset)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::facet_wrap(~ plot_type, scales = "free") + 
      ggplot2::labs(title = paste("Overall Distribution"), 
         x = paste("Log Gene Expression", gene_name), 
         y = "Density") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                     plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10))
  
  return(plot_obj)
}