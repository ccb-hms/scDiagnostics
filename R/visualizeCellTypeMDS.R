#' Visualizing Reference and Query Cell Types using MDS
#'
#' This function facilitates the assessment of similarity between
#' reference and query datasets through Multidimensional Scaling (MDS)
#' scatter plots. It allows the visualization of cell types,
#' color-coded with user-defined custom colors, based on a
#' dissimilarity matrix computed from a user-selected gene set.
#' 
#' @details To evaluate dataset similarity, the function selects
#'     specific subsets of cells from both reference and query
#'     datasets. It then calculates Spearman correlations between gene
#'     expression profiles, deriving a dissimilarity matrix. This
#'     matrix undergoes Classical Multidimensional Scaling (MDS) for
#'     visualization, presenting cell types in a scatter plot,
#'     distinguished by colors defined by the user.
#' 
#' @param query_data A \code{\linkS4class{SingleCellExperiment}}
#'     containing the single-cell expression data and metadata.
#'
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}}
#'     object containing the single-cell expression data and metadata.
#'
#' @param mdata A character vector representing cell types from both
#'     the query and reference datasets.  These cell types are labels
#'     used for visualization in the Multidimensional Scaling (MDS)
#'     plot.  The vector should include names of cell types as they
#'     appear in both datasets, prefixed appropriately to distinguish
#'     between query and reference (e.g., "Query CD4",
#'     "Reference CD8").
#' 
#' @param cell_type_colors A named vector of colors corresponding to
#'     the cell types specified in mdata.  Each color is used to
#'     uniquely identify a cell type in the MDS plot.
#' 
#' @param legend_order A character vector specifying the order of cell
#'     types in the plot legend.
#'
#' @return A ggplot object representing the MDS scatter plot with cell
#'     type coloring.
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
#' cell_types <- c("Query CD8", 
#'                 "Reference CD8", 
#'                  "Query CD4", 
#'                  "Reference CD4", 
#'                  "Query B_and_plasma", 
#'                  "Reference B_and_plasma")
#'                  
#' legend_order <- cell_types
#'
#' # Define the colors for cell types
#' color_palette <- brewer.pal(length(cell_types), "Paired")
#' color_mapping <- setNames(color_palette, cell_types)
#' cell_type_colors <- color_mapping[cell_types]
#'
#' # Generate the MDS scatter plot with cell type coloring
#' plot <- visualizeCellTypeMDS(query_data = query_data_subset, 
#'                             reference_data = ref_data_subset, 
#'                             mdata = mdata, 
#'                             cell_type_colors = cell_type_colors, 
#'                             legend_order = legend_order)
#' print(plot)
#'
#' @importFrom stats cmdscale cor
#' @importFrom ggplot2 ggplot
#' @importFrom SummarizedExperiment assay
#'
#' @export
#' 
visualizeCellTypeMDS <- function(query_data, 
                                 reference_data, 
                                 mdata, 
                                 cell_type_colors, 
                                 legend_order) {

    ## Sanity checks
  
    ## Check if query_data is a SingleCellExperiment object
    if (!is(query_data, "SingleCellExperiment")) {
        stop("query_data must be a SingleCellExperiment object.")
    }
  
    ## Check if reference_data is a SingleCellExperiment object
    if (!is(reference_data, "SingleCellExperiment")) {
        stop("reference_data must be a SingleCellExperiment object.")
    }
  
    ## Extract logcounts
    queryExp <- as.matrix(assay(query_data, "logcounts"))
    refExp <- as.matrix(assay(reference_data, "logcounts"))

    ## Compute correlation and dissimilarity matrix
    df <- cbind(queryExp, refExp)
    corMat <- cor(df, method = "spearman")
    disMat <- (1 - corMat)
    cmd <- cmdscale(disMat)

    ## Create the plot object
    matx <- data.frame(
        Dim1 = cmd[, 1],
        Dim2 = cmd[, 2],
        Type = factor(mdata, levels = legend_order)
    )
    
    plot <- ggplot(matx, aes(x = Dim1, y = Dim2, color = Type)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = cell_type_colors) +
        theme_bw() +
        guides(color = guide_legend(title = "Cell Types"))
    
    plot
}
