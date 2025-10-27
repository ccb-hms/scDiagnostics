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
#' Multiple normalization options are available:
#' - "z_score": Standard z-score normalization within each dataset
#' - "min_max": Min-max scaling to [0,1] range within each dataset
#' - "rank": Maps values to quantile ranks (0-100 scale)
#' - "none": No transformation (preserves original scale differences)
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_type A cell type to plot (e.g., c("T-cell", "B-cell")).
#' @param gene_name The gene name for which the distribution is to be visualized.
#' @param normalization Method for normalizing expression values. Options: "z_score" (default), "min_max", "rank", "none".
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells_query Maximum number of query cells to retain after cell type filtering. If NULL,
#' no downsampling of query cells is performed. Default is NULL.
#' @param max_cells_ref Maximum number of reference cells to retain after cell type filtering. If NULL,
#' no downsampling of reference cells is performed. Default is NULL.
#'
#' @return A ggplot object containing density plots comparing reference and query distributions.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
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
#'                      query_cell_type_col = c("expert_annotation", "SingleR_annotation")[1],
#'                      gene_name = "CD8A",
#'                      cell_type = "CD4",
#'                      normalization = "z_score")
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom stats sd
#' @import SingleCellExperiment
#'
# Function to plot the expression of a marker
plotMarkerExpression <- function(query_data,
                                 reference_data,
                                 ref_cell_type_col,
                                 query_cell_type_col,
                                 cell_type,
                                 gene_name,
                                 assay_name = "logcounts",
                                 normalization = c("z_score", "min_max", "rank", "none"),
                                 max_cells_query = NULL,
                                 max_cells_ref = NULL) {

    # Match normalization argument
    normalization <- match.arg(normalization)

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  assay_name = assay_name,
                  max_cells_query = max_cells_query,
                  max_cells_ref = max_cells_ref)

    # Convert cell type columns to character if needed
    query_data <- convertColumnsToCharacter(sce_object = query_data,
                                            convert_cols = query_cell_type_col)
    reference_data <- convertColumnsToCharacter(sce_object = reference_data,
                                                convert_cols = ref_cell_type_col)

    # Select cell types
    cell_type <- selectCellTypes(query_data = query_data,
                                 reference_data = reference_data,
                                 query_cell_type_col = query_cell_type_col,
                                 ref_cell_type_col = ref_cell_type_col,
                                 cell_types = cell_type,
                                 dual_only = TRUE,
                                 n_cell_types = NULL)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce_object = query_data,
                                cell_type_col = query_cell_type_col,
                                max_cells = max_cells_query)
    reference_data <- downsampleSCE(sce_object = reference_data,
                                    cell_type_col = ref_cell_type_col,
                                    max_cells = max_cells_ref)

    # Check if gene_name is present in both query_data and reference_data
    if (!(gene_name %in% rownames(assay(query_data)) &&
          gene_name %in% rownames(assay(reference_data)))) {
        stop("gene_name: \'",
             gene_name,
             "\' is not present in the row names of both \'query_data\' and \'reference_data\'.")
    }

    # Get expression of the specified gene for reference and query datasets
    ref_gene_expression <- assay(reference_data, assay_name)[gene_name, ]
    query_gene_expression <- assay(query_data, assay_name)[gene_name, ]
    ref_gene_expression_specific <- assay(
        reference_data, assay_name)[gene_name,
                                    which(reference_data[[ref_cell_type_col]] %in%
                                              cell_type)]
    query_gene_expression_specific <- assay(
        query_data, assay_name)[gene_name,
                                which(query_data[[query_cell_type_col]] %in%
                                          cell_type)]

    # Transformation functions
    .quantileTransformation <- function(x) {
        return(100 * rank(x, ties.method = "average") / length(x))
    }

    .zScoreTransformation <- function(x) {
        if(sd(x) == 0) return(rep(0, length(x)))
        return((x - mean(x)) / sd(x))
    }

    .minMaxTransformation <- function(x) {
        if(max(x) == min(x)) return(rep(0.5, length(x)))
        return((x - min(x)) / (max(x) - min(x)))
    }

    # Apply selected normalization
    if (normalization == "rank") {
        ref_gene_expression_norm <-
            .quantileTransformation(ref_gene_expression)
        query_gene_expression_norm <-
            .quantileTransformation(query_gene_expression)
        ref_gene_expression_specific_norm <-
            .quantileTransformation(ref_gene_expression_specific)
        query_gene_expression_specific_norm <-
            .quantileTransformation(query_gene_expression_specific)
        x_label <-
            paste("Quantile Rank Normalized Gene Expression:", gene_name)

    } else if (normalization == "z_score") {
        ref_gene_expression_norm <-
            .zScoreTransformation(ref_gene_expression)
        query_gene_expression_norm <-
            .zScoreTransformation(query_gene_expression)
        ref_gene_expression_specific_norm <-
            .zScoreTransformation(ref_gene_expression_specific)
        query_gene_expression_specific_norm <-
            .zScoreTransformation(query_gene_expression_specific)
        x_label <-
            paste("Z-Score Normalized Gene Expression:", gene_name)

    } else if (normalization == "min_max") {
        ref_gene_expression_norm <-
            .minMaxTransformation(ref_gene_expression)
        query_gene_expression_norm <-
            .minMaxTransformation(query_gene_expression)
        ref_gene_expression_specific_norm <-
            .minMaxTransformation(ref_gene_expression_specific)
        query_gene_expression_specific_norm <-
            .minMaxTransformation(query_gene_expression_specific)
        x_label <- paste("Min-Max Normalized Gene Expression:", gene_name)
    } else if (normalization == "none") {
        ref_gene_expression_norm <-
            ref_gene_expression
        query_gene_expression_norm <-
            query_gene_expression
        ref_gene_expression_specific_norm <-
            ref_gene_expression_specific
        query_gene_expression_specific_norm <-
            query_gene_expression_specific
        x_label <-
            paste("Log-Normalized Gene Expression:", gene_name)

    }

    # Create a combined vector of gene expression values
    combined_gene_expression <- c(
        ref_gene_expression_norm,
        query_gene_expression_norm,
        ref_gene_expression_specific_norm,
        query_gene_expression_specific_norm)

    # Create a grouping vector for dataset types
    dataset_types <- rep(c("Reference", "Query", "Reference", "Query"),
                         times = c(
                             length(ref_gene_expression),
                             length(query_gene_expression),
                             length(ref_gene_expression_specific),
                             length(query_gene_expression_specific)))

    # Combine the gene expression values and dataset types
    marker_data <- data.frame(
        GeneExpression = combined_gene_expression,
        Dataset = dataset_types,
        plot_type = rep(c("Overall Distribution", "Cell Type-Specific Distribution"),
                        times = c(length(ref_gene_expression) +
                                      length(query_gene_expression),
                                  length(ref_gene_expression_specific) +
                                      length(query_gene_expression_specific))))
    marker_data[["Dataset"]] <- factor(marker_data[["Dataset"]],
                                       levels = c("Query", "Reference"))


    # Create a stacked density plot
    plot_obj <- ggplot2::ggplot(
        marker_data,
        ggplot2::aes(x = .data[["GeneExpression"]],
                     y = .data[["Dataset"]],
                     fill = .data[["Dataset"]])) +
        ggridges::geom_density_ridges(
            alpha = 0.7,
            scale = 1,
            rel_min_height = 0.01
        ) +
        ggplot2::scale_y_discrete(
            expand = ggplot2::expansion(mult = c(0.02, 0.02))
        ) +
        ggplot2::scale_fill_manual(
            values = c("Query" = "#B565D8", "Reference" = "#5A9BD8"),
            name = "Dataset"
        ) +
        ggplot2::facet_wrap(~ .data[["plot_type"]], scales = "free") +
        ggplot2::labs(title = NULL,
                      x = x_label,
                      y = "") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            panel.grid.major.y = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "white",
                                                     color = "black"),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank())
    return(plot_obj)
}
