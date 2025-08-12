#' @title Calculate Top Loading Gene Expression Shifts
#'
#' @description
#' This function identifies genes with the highest loadings for specified principal components
#' and performs statistical tests to detect distributional differences between query and reference data.
#' It also calculates the proportion of variance explained by each principal component within
#' specific cell types.
#'
#' @details
#' This function extracts the top loading genes for each specified principal component from the reference
#' PCA space and performs distributional comparisons between query and reference data. For each gene,
#' it performs statistical tests to identify genes that may be causing PC-specific alignment issues
#' between datasets. A key feature is the calculation of cell-type-specific variance explained by
#' global PCs, providing a more nuanced view of how major biological axes affect individual populations.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to analyze. If NULL, all common cell types are used.
#' @param pc_subset A numeric vector specifying which principal components to analyze. Default is 1:5.
#' @param n_top_loadings Number of top loading genes to analyze per PC. Default is 50.
#' @param p_value_threshold P-value threshold for statistical significance. Default is 0.05.
#' @param adjust_method Method for multiple testing correction. Default is "fdr".
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A list containing:
#' \itemize{
#'   \item PC results: Named elements for each PC (e.g., "PC1", "PC2") containing data frames with gene-level analysis results.
#'   \item expression_data: Matrix of expression values for all analyzed genes (genes × cells).
#'   \item cell_metadata: Data frame with columns: cell_id, dataset, cell_type, original_index.
#'   \item gene_metadata: Data frame with columns: gene, pc, loading for all analyzed genes.
#'   \item percent_var: Named numeric vector of global percent variance explained for each analyzed PC.
#'   \item cell_type_variance: A data frame detailing the percent of variance a global PC explains within specific cell types for both query and reference datasets.
#' }
#'
#' The `cell_type_variance` data frame contains columns: pc, cell_type, dataset, percent_variance.
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateTopLoadingGeneShiftsObject}}
#'
#' @importFrom stats wilcox.test var p.adjust na.omit setNames
#'
calculateTopLoadingGeneShifts <- function(query_data,
                                          reference_data,
                                          query_cell_type_col,
                                          ref_cell_type_col,
                                          cell_types = NULL,
                                          pc_subset = 1:5,
                                          n_top_loadings = 50,
                                          p_value_threshold = 0.05,
                                          adjust_method = "fdr",
                                          assay_name = "logcounts",
                                          max_cells = 2500) {

    # Input validation
    if (!is.numeric(n_top_loadings) || length(n_top_loadings) != 1 || n_top_loadings <= 0) {
        stop("n_top_loadings must be a positive integer")
    }

    if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1 ||
        p_value_threshold < 0 || p_value_threshold > 1) {
        stop("p_value_threshold must be between 0 and 1")
    }

    # Check if cell type columns exist
    if (!query_cell_type_col %in% names(SummarizedExperiment::colData(query_data))) {
        stop(paste("Column '", query_cell_type_col, "' not found in query_data colData"))
    }

    if (!ref_cell_type_col %in% names(SummarizedExperiment::colData(reference_data))) {
        stop(paste("Column '", ref_cell_type_col, "' not found in reference_data colData"))
    }

    # Get common cell types
    if (is.null(cell_types)) {
        cell_types <- stats::na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                              query_data[[query_cell_type_col]])))
    }

    # Get PCA results
    pca_attrs <- attributes(SingleCellExperiment::reducedDim(reference_data, "PCA"))
    pca_rotation <- pca_attrs[["rotation"]]
    percent_var <- pca_attrs[["percentVar"]]
    if (is.null(pca_rotation) || is.null(percent_var)) {
        stop("PCA rotation matrix and/or percentVar not found in reference_data's reducedDim 'PCA' attributes.")
    }
    if (max(pc_subset) > ncol(pca_rotation)) {
        stop("pc_subset contains indices greater than the number of PCs available in the reference.")
    }
    pc_percent_var <- percent_var[pc_subset]
    names(pc_percent_var) <- paste0("PC", pc_subset)

    # Get common genes
    common_genes <- intersect(rownames(query_data), rownames(reference_data))
    common_genes <- intersect(common_genes, rownames(pca_rotation))
    if (length(common_genes) == 0) {
        stop("No common genes found between query, reference, and PCA rotation matrix.")
    }
    pca_rotation <- pca_rotation[common_genes, ]

    # Prepare cell indices
    query_cell_indices <- split(seq_len(ncol(query_data)), query_data[[query_cell_type_col]])
    ref_cell_indices <- split(seq_len(ncol(reference_data)), reference_data[[ref_cell_type_col]])
    query_cell_indices <- query_cell_indices[names(query_cell_indices) %in% cell_types]
    ref_cell_indices <- ref_cell_indices[names(ref_cell_indices) %in% cell_types]
    available_cell_types <- intersect(names(query_cell_indices), names(ref_cell_indices))
    if (length(available_cell_types) == 0) {
        warning("No common cell types with sufficient cells found between datasets.")
        return(list())
    }

    # Gene-level statistical analysis
    all_top_genes <- character(0)
    gene_metadata_list <- list()
    pc_results <- stats::setNames(vector("list", length(pc_subset)), paste0("PC", pc_subset))

    for (pc in pc_subset) {
        pc_name <- paste0("PC", pc)
        pc_loadings <- pca_rotation[, pc]
        top_loading_indices <- order(abs(pc_loadings), decreasing = TRUE)[1:min(n_top_loadings, length(pc_loadings))]
        top_genes <- names(pc_loadings)[top_loading_indices]
        top_loadings_vals <- pc_loadings[top_loading_indices]

        all_top_genes <- unique(c(all_top_genes, top_genes))
        gene_metadata_list[[pc_name]] <- data.frame(gene = top_genes, pc = pc, loading = top_loadings_vals, stringsAsFactors = FALSE)

        pc_result_list <- list()
        for (ct in available_cell_types) {
            query_cells_ct <- query_cell_indices[[ct]]
            ref_cells_ct <- ref_cell_indices[[ct]]
            if (length(query_cells_ct) < 3 || length(ref_cells_ct) < 3) next

            query_expr_ct <- SummarizedExperiment::assay(query_data, assay_name)[top_genes, query_cells_ct, drop = FALSE]
            ref_expr_ct <- SummarizedExperiment::assay(reference_data, assay_name)[top_genes, ref_cells_ct, drop = FALSE]

            gene_results <- processGenesSimple(top_genes, top_loadings_vals, query_expr_ct, ref_expr_ct, ct)
            if (length(gene_results) > 0) pc_result_list <- c(pc_result_list, gene_results)
        }

        if (length(pc_result_list) > 0) {
            df <- do.call(rbind, pc_result_list)
            df$p_adjusted <- stats::p.adjust(df$p_value, method = adjust_method)
            df$significant <- df$p_adjusted <= p_value_threshold
            pc_results[[pc_name]] <- df[order(df$p_adjusted), ]
            rownames(pc_results[[pc_name]]) <- NULL
        } else {
            pc_results[[pc_name]] <- data.frame()
        }
    }
    gene_metadata <- do.call(rbind, gene_metadata_list)
    rownames(gene_metadata) <- NULL

    # Cell-type-specific variance explained by global PCs
    var_explained_list <- list()
    full_ref_assay <- SummarizedExperiment::assay(reference_data, assay_name)[common_genes, ]
    full_query_assay <- SummarizedExperiment::assay(query_data, assay_name)[common_genes, ]

    for (ct in available_cell_types) {
        ref_cells_ct <- ref_cell_indices[[ct]]
        query_cells_ct <- query_cell_indices[[ct]]

        ref_expr_ct <- full_ref_assay[, ref_cells_ct, drop = FALSE]
        query_expr_ct <- full_query_assay[, query_cells_ct, drop = FALSE]

        # CHANGE: Replaced matrixStats::rowVars with base R apply()
        total_var_ref <- if(ncol(ref_expr_ct) > 1) sum(apply(as.matrix(ref_expr_ct), 1, stats::var)) else 0
        total_var_query <- if(ncol(query_expr_ct) > 1) sum(apply(as.matrix(query_expr_ct), 1, stats::var)) else 0

        for (pc in pc_subset) {
            pc_name <- paste0("PC", pc)
            loadings_pc <- pca_rotation[, pc, drop = FALSE]

            if (total_var_ref > .Machine$double.eps) {
                ref_scores <- crossprod(ref_expr_ct, loadings_pc)
                var_scores_ref <- stats::var(as.vector(ref_scores))
                pct_var_ref <- (var_scores_ref / total_var_ref) * 100
                var_explained_list[[length(var_explained_list) + 1]] <- data.frame(
                    pc = pc_name, cell_type = ct, dataset = "Reference", percent_variance = pct_var_ref, stringsAsFactors = FALSE)
            }
            if (total_var_query > .Machine$double.eps) {
                query_scores <- crossprod(query_expr_ct, loadings_pc)
                var_scores_query <- stats::var(as.vector(query_scores))
                pct_var_query <- (var_scores_query / total_var_query) * 100
                var_explained_list[[length(var_explained_list) + 1]] <- data.frame(
                    pc = pc_name, cell_type = ct, dataset = "Query", percent_variance = pct_var_query, stringsAsFactors = FALSE)
            }
        }
    }
    cell_type_variance_df <- do.call(rbind, var_explained_list)
    rownames(cell_type_variance_df) <- NULL

    # Data for plotting
    query_relevant_cells <- unlist(query_cell_indices[available_cell_types], use.names = FALSE)
    ref_relevant_cells <- unlist(ref_cell_indices[available_cell_types], use.names = FALSE)

    query_expr_plot <- SummarizedExperiment::assay(query_data, assay_name)[all_top_genes, query_relevant_cells, drop = FALSE]
    ref_expr_plot <- SummarizedExperiment::assay(reference_data, assay_name)[all_top_genes, ref_relevant_cells, drop = FALSE]
    expression_data <- cbind(ref_expr_plot, query_expr_plot)

    cell_metadata <- rbind(
        data.frame(cell_id = colnames(ref_expr_plot), dataset = "Reference",
                   cell_type = reference_data[[ref_cell_type_col]][ref_relevant_cells],
                   original_index = ref_relevant_cells, stringsAsFactors = FALSE),
        data.frame(cell_id = colnames(query_expr_plot), dataset = "Query",
                   cell_type = query_data[[query_cell_type_col]][query_relevant_cells],
                   original_index = query_relevant_cells, stringsAsFactors = FALSE)
    )

    # Assemble final results
    final_results <- c(
        pc_results,
        list(
            expression_data = expression_data,
            cell_metadata = cell_metadata,
            gene_metadata = gene_metadata,
            pca_rotation = pca_rotation,
            percent_var = pc_percent_var,
            cell_type_variance = cell_type_variance_df
        )
    )

    class(final_results) <- c(class(final_results), "calculateTopLoadingGeneShiftsObject")
    return(final_results)
}

#' @title Process Genes for Statistical Analysis
#'
#' @description
#' Internal helper function that performs statistical testing on top loading genes
#' to compare their expression distributions between query and reference datasets.
#'
#' @param top_genes Character vector of gene names to analyze.
#' @param top_loadings Numeric vector of PC loading values corresponding to the genes.
#' @param query_expr_matrix Numeric matrix of expression values for query data (genes × cells).
#' @param ref_expr_matrix Numeric matrix of expression values for reference data (genes × cells).
#' @param cell_type Character string specifying the cell type being analyzed.
#'
#' @keywords internal
#' @return A list of data frames, where each data frame contains results for one gene.
#'
processGenesSimple <- function(top_genes, top_loadings,
                               query_expr_matrix, ref_expr_matrix,
                               cell_type) {

    result_list <- list()

    for (i in seq_along(top_genes)) {
        gene <- top_genes[i]
        loading <- top_loadings[i]

        query_expr <- query_expr_matrix[gene, ]
        ref_expr <- ref_expr_matrix[gene, ]

        query_expr <- query_expr[!is.na(query_expr)]
        ref_expr <- ref_expr[!is.na(ref_expr)]

        if (length(query_expr) < 3 || length(ref_expr) < 3) {
            next
        }

        mean_query <- mean(query_expr)
        mean_reference <- mean(ref_expr)

        test_result <- suppressWarnings(
            stats::wilcox.test(query_expr, ref_expr, alternative = "two.sided")
        )

        result_row <- data.frame(
            gene = gene,
            loading = loading,
            cell_type = cell_type,
            p_value = test_result[["p.value"]],
            mean_query = mean_query,
            mean_reference = mean_reference,
            stringsAsFactors = FALSE
        )

        result_list[[length(result_list) + 1]] <- result_row
    }

    return(result_list)
}
