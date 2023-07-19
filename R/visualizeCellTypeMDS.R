#' Visualizing Reference and Query Cell Types using MDS
#'
#' This function performs Multidimensional Scaling (MDS) on the query and reference datasets and generates a scatter plot
#' with cell types color-coded using custom colors provided by the user. The distance matrix used for MDS is computed
#' based on a user choice of genes and cells.
#'
#' @param query_data An object of class "SingleCellExperiment" containing a numeric expression matrix.
#' @param ref_data An object of class "SingleCellExperiment" containing a numeric expression matrix.
#' @param mdata A vector of cell types from the query and reference datasets to be visualized in low-dimensional space.
#' @param colors A vector of custom colors to be used for the cell types.
#' @param legend_order A vector specifying the desired order of the legend items.
#'
#' @importFrom stats cmdscale cor
#' @importFrom ggplot2 ggplot
#'
#' @return A ggplot object representing the MDS scatter plot with cell type coloring.
#' @export
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(RColorBrewer)
#'
#' # load data
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
#' # Selecting highly variable genes
#' ref_var <- getTopHVGs(ref_data, n=2000)
#' query_var <- getTopHVGs(query_data, n=2000)
#'
#' # Intersect the gene symbols to obtain common genes
#' common_genes <- intersect(ref_var, query_var)
#'
#' # Select desired cell types
#' selected_cell_types <- c("CD4", "CD8", "B_and_plasma")
#' ref_data_subset <- ref_data[common_genes, ref_data$reclustered.broad %in% selected_cell_types]
#' query_data_subset <- query_data[common_genes, query_data$reclustered.broad %in% selected_cell_types]
#'
#' # Extract cell types for visualization
#' ref_labels <- ref_data_subset$reclustered.broad
#' query_labels <- query_data_subset$reclustered.broad
#'
#' # Combine the cell type labels from both datasets
#' mdata <- c(paste("Query", query_labels), paste("Reference", ref_labels))
#'
#' # Define the cell types and legend order
#' cell_types <- c("Query CD8", "Reference CD8", "Query CD4", "Reference CD4", "Query B_and_plasma", "Reference B_and_plasma")
#' legend_order <- cell_types
#'
#' # Define the colors for cell types
#' color_palette <- brewer.pal(length(cell_types), "Paired")
#' color_mapping <- setNames(color_palette, cell_types)
#' cell_type_colors <- color_mapping[cell_types]
#'
#' # Generate the MDS scatter plot with cell type coloring
#' plot <- visualizeCellTypeMDS(query_data_subset, ref_data_subset, mdata, cell_type_colors, legend_order)
#' print(plot)
#'
visualizeCellTypeMDS <- function(query_data, ref_data, mdata, colors, legend_order) {

  queryExp <- as.matrix(assay(query_data, "logcounts"))
  refExp <- as.matrix(assay(ref_data, "logcounts"))

  # Compute correlation and dissimilarity matrix
  df <- cbind(queryExp, refExp)
  corMat <- cor(df, method = "spearman")
  disMat <- (1 - corMat)
  cmd <- cmdscale(disMat)

  # Create the plot object
  matx <- data.frame(Dim1 = cmd[, 1], Dim2 = cmd[, 2], Type = factor(mdata, levels = legend_order))

  plot <- ggplot(matx, aes(x = Dim1, y = Dim2, color = Type)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = colors) +
    theme_bw() +
    guides(color = guide_legend(title = "Cell Types"))

  return(plot)
}
