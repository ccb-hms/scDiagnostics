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
#' To make the gene expression scales comparable between the datasets, the gene expression values 
#' are transformed using z-rank normalization. This transformation ranks the expression values 
#' and then scales the ranks to have a mean of 0 and a standard deviation of 1, which helps 
#' in standardizing the distributions for comparison.
#' 
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_type A vector of cell type cell_types to plot (e.g., c("T-cell", "B-cell")).
#' @param gene_name The gene name for which the distribution is to be visualized.
#'
#' @return A gtable object containing two arranged density plots as grobs. 
#'         The first plot shows the overall gene expression distribution, 
#'         and the second plot displays the cell type-specific expression 
#'         distribution.
#'         
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#' 
#' # Note: Users can use SingleR or any other method to obtain the cell type annotations.
#' plotMarkerExpression(reference_data = reference_data, 
#'                      query_data = query_data, 
#'                      ref_cell_type_col = "expert_annotation", 
#'                      query_cell_type_col = "SingleR_annotation", 
#'                      gene_name = "VPREB3", 
#'                      cell_type = "B_and_plasma")
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom stats sd
#' @import SingleCellExperiment
#' @export
plotMarkerExpression <- function(reference_data, 
                                 query_data, 
                                 ref_cell_type_col, 
                                 query_cell_type_col, 
                                 cell_type,
                                 gene_name) {
    
    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_type)
    
    # Get common cell types if they are not specified by user
    if(is.null(cell_type)){
        cell_type <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }
    
    # Check if gene_name is present in both query_data and reference_data
    if (!(gene_name %in% rownames(assay(query_data)) && 
          gene_name %in% rownames(assay(reference_data)))) {
        stop("gene_name: '", gene_name, "' is not present in the row names of both \'query_data\' and \'reference_data\'.")
    }
    
    # Get expression of the specified gene for reference and query datasets
    ref_gene_expression <- assay(reference_data, "logcounts")[gene_name, ]
    query_gene_expression <- assay(query_data, "logcounts")[gene_name, ]
    ref_gene_expression_specific <- assay(reference_data, "logcounts")[gene_name, which(reference_data[[ref_cell_type_col]] %in% cell_type)]
    query_gene_expression_specific <- assay(query_data, "logcounts")[gene_name, which(query_data[[query_cell_type_col]] %in% cell_type)]
    
    # Z-rank transformation
    .rankTransformation <- function(x) {
        ranks <- rank(x, ties.method = "average")
        z_ranks <- (ranks - mean(ranks)) / sd(ranks)
        return(z_ranks)
    }
    ref_gene_expression_zr <- .rankTransformation(ref_gene_expression)
    query_gene_expression_zr <- .rankTransformation(query_gene_expression)
    ref_gene_expression_specific_zr <- .rankTransformation(ref_gene_expression_specific)
    query_gene_expression_specific_zr <- .rankTransformation(query_gene_expression_specific)
    
    # Create a combined vector of gene expression values
    combined_gene_expression <- c(ref_gene_expression_zr, query_gene_expression_zr,
                                  ref_gene_expression_specific_zr, query_gene_expression_specific_zr)
    
    # Create a grouping vector for dataset cell_types
    dataset_cell_types <- rep(c("Reference", "Query", "Reference", "Query"), 
                               times = c(length(ref_gene_expression), length(query_gene_expression),
                                         length(ref_gene_expression_specific), length(query_gene_expression_specific)))
    
    # Combine the gene expression values and dataset cell_types
    data <- data.frame(GeneExpression = combined_gene_expression,
                       Dataset = dataset_cell_types,
                       plot_type = rep(c("Overall Distribution", "Cell Type-Specific Distribution"), 
                                       times = c(length(ref_gene_expression) + length(query_gene_expression),
                                                 length(ref_gene_expression_specific) + length(query_gene_expression_specific))))
    
    # Create a stacked density plot using ggplot2 for overall dataset
    plot_obj <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["GeneExpression"]], fill = .data[["Dataset"]])) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::facet_wrap(~ .data[["plot_type"]], scales = "free") + 
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
