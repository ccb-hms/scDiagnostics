#' Visualizing Reference and Query Cell Types using MDS
#'
#' This function facilitates the assessment of similarity between reference and query datasets 
#' through Multidimensional Scaling (MDS) scatter plots. It allows the visualization of cell types, 
#' color-coded with user-defined custom colors, based on a dissimilarity matrix computed from a 
#' user-selected gene set.
#' 
#' @details To evaluate dataset similarity, the function selects specific subsets of cells from 
#' both reference and query datasets. It then calculates Spearman correlations between gene expression profiles, 
#' deriving a dissimilarity matrix. This matrix undergoes Classical Multidimensional Scaling (MDS) for 
#' visualization, presenting cell types in a scatter plot, distinguished by colors defined by the user.
#' 
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} containing the single-cell 
#' expression data and metadata.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing the single-cell 
#' expression data and metadata.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data} 
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data} 
#' that identifies the cell types.
#' 
#' @return A ggplot object representing the MDS scatter plot with cell type coloring.
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
#' set.seed(100)
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
#' ref_data_subset <- ref_data[common_genes,]
#' query_data_subset <- query_data[common_genes,]
#' 
#' # Generate the MDS scatter plot with cell type coloring
#' plot <- visualizeCellTypeMDS(query_data = query_data_subset, 
#'                              reference_data = ref_data_subset, 
#'                              cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid")[1:4],
#'                              query_cell_type_col = "labels", 
#'                              ref_cell_type_col = "reclustered.broad")
#' print(plot)
#'
#' @importFrom stats cmdscale cor
#' @importFrom ggplot2 ggplot
#' @importFrom SummarizedExperiment assay
#' @export
#' 
visualizeCellTypeMDS <- function(query_data, 
                                 reference_data, 
                                 cell_types = NULL,
                                 query_cell_type_col = NULL, 
                                 ref_cell_type_col = NULL) {

    # Check if query_data is a SingleCellExperiment object
    if (!is(query_data, "SingleCellExperiment")) {
    stop("query_data must be a SingleCellExperiment object.")
    }
    
    # Check if reference_data is a SingleCellExperiment object
    if (!is(reference_data, "SingleCellExperiment")) {
    stop("reference_data must be a SingleCellExperiment object.")
    }
    
    # Check if query_cell_type_col is a valid column name in query_data
    if (!query_cell_type_col %in% colnames(colData(query_data))) {
      stop("query_cell_type_col: '", query_cell_type_col, "' is not a valid column name in query_data.")
    }
    
    # Check if ref_cell_type_col is a valid column name in reference_data
    if (!ref_cell_type_col %in% colnames(colData(reference_data))) {
      stop("ref_cell_type_col: '", ref_cell_type_col, "' is not a valid column name in reference_data.")
    }
    
    # Check if cell types available in both single-cell experiments
    if(!all(cell_types %in% reference_data[[ref_cell_type_col]]) || 
       !all(cell_types %in% query_data[[query_cell_type_col]]))
        stop("One or more of the specified cell types are not available in \'reference_data\' or \'query_data\'.")
    
    # Subset data
    query_data <- query_data[, which(query_data[[query_cell_type_col]] %in% cell_types)]
    reference_data <- reference_data[, which(reference_data[[ref_cell_type_col]] %in% cell_types)]
    
    # Extract logcounts
    queryExp <- as.matrix(assay(query_data, "logcounts"))
    refExp <- as.matrix(assay(reference_data, "logcounts"))
    
    # Compute correlation and dissimilarity matrix
    df <- cbind(queryExp, refExp)
    corMat <- cor(df, method = "spearman")
    disMat <- (1 - corMat)
    cmd <- data.frame(cmdscale(disMat), c(rep("Query", ncol(queryExp)), rep("Reference", ncol(refExp))),
                      c(query_data[[query_cell_type_col]], reference_data[[ref_cell_type_col]]))
    colnames(cmd) <- c("Dim1", "Dim2", "dataset", "cellType")
    cmd <- na.omit(cmd)
    cmd$cell_type_dataset <- paste(cmd$dataset, cmd$cellType, sep = " ")

    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)), rep(sort(cell_types), each = 2))
    cmd$cell_type_dataset <- factor(cmd$cell_type_dataset, levels = order_combinations)
    
    # Define the colors for cell types
    color_mapping <- setNames(RColorBrewer::brewer.pal(length(order_combinations), "Paired"), order_combinations)
    cell_type_colors <- color_mapping[order_combinations]
    
    plot <- ggplot2::ggplot(cmd, aes(x = Dim1, y = Dim2, color = cell_type_dataset)) +
      ggplot2::geom_point(alpha = 0.5, size = 1) +
      ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Types") + 
      ggplot2::theme_bw() +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Cell Types"))
    return(plot)
}
