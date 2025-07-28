#' @title Calculate Top Loading Gene Expression Shifts
#'
#' @description
#' This function identifies genes with the highest loadings for specified principal components
#' and performs statistical tests to detect distributional differences between query and reference data.
#'
#' @details
#' This function extracts the top loading genes for each specified principal component from the reference
#' PCA space and performs distributional comparisons between query and reference data. For each gene,
#' it performs statistical tests to identify genes that may be causing PC-specific alignment issues
#' between datasets.
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
#'   \item PC results: Named elements for each PC (e.g., "PC1", "PC2") containing data frames with analysis results
#'   \item expression_data: Matrix of expression values for all analyzed genes (genes × cells)
#'   \item cell_metadata: Data frame with columns: cell_id, dataset, cell_type, original_index
#'   \item gene_metadata: Data frame with columns: gene, pc, loading for all analyzed genes
#'   \item percent_var: Named numeric vector of percent variance explained for each analyzed PC
#' }
#'
#' Each PC data frame contains columns:
#' \itemize{
#'   \item gene: Gene symbol
#'   \item loading: PC loading value for the gene
#'   \item cell_type: Cell type analyzed
#'   \item p_value: Raw p-value from Wilcoxon rank-sum test
#'   \item p_adjusted: Adjusted p-value
#'   \item mean_query: Mean expression in query data
#'   \item mean_reference: Mean expression in reference data
#'   \item significant: Logical indicating statistical significance
#' }
#'
#' @export
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @seealso \code{\link{plot.calculateTopLoadingGeneShiftsObject}}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Compute distributional shifts for genes with top loadings
#' gene_shifts <- calculateTopLoadingGeneShifts(reference_data = reference_data,
#'                                              query_data = query_data,
#'                                              query_cell_type_col = "SingleR_annotation",
#'                                              ref_cell_type_col = "expert_annotation",
#'                                              cell_types = NULL,
#'                                              pc_subset = 1:3,
#'                                              n_top_loadings = 50,
#'                                              assay_name = "logcounts",
#'                                              p_value_threshold = 0.05,
#'                                              adjust_method = "fdr")
#'
#' # Plot gene shifts
#' plot(gene_shifts,
#'      cell_type = "CD8",
#'      pc_subset = 1:3,
#'      plot_by = "p_adjusted",
#'      n_genes = 10,
#'      significance_threshold = 0.05)
#'
#' @importFrom stats wilcox.test
#'
# Function to calculate distributional shifts for genes with top loadings
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

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)

    # Additional input validation
    if (!is.numeric(n_top_loadings) || n_top_loadings <= 0) {
        stop("n_top_loadings must be a positive integer")
    }
    if (!is.numeric(p_value_threshold) || p_value_threshold <= 0 || p_value_threshold > 1) {
        stop("p_value_threshold must be between 0 and 1")
    }

    # Get common cell types if not specified
    if (is.null(cell_types)) {
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Get PCA rotation matrix (loadings) and percent variance
    pca_rotation <- attributes(reducedDim(reference_data, "PCA"))[["rotation"]]
    percent_var <- attributes(reducedDim(reference_data, "PCA"))[["percentVar"]]

    if (is.null(pca_rotation)) {
        stop("PCA rotation matrix not found in reference data")
    }

    if (is.null(percent_var)) {
        stop("Percent variance not found in reference data PCA")
    }

    # Extract percent variance for analyzed PCs
    pc_percent_var <- percent_var[pc_subset]
    names(pc_percent_var) <- paste0("PC", pc_subset)

    # Get common genes between datasets
    common_genes <- intersect(rownames(query_data), rownames(reference_data))
    common_genes <- intersect(common_genes, rownames(pca_rotation))

    if (length(common_genes) == 0) {
        stop("No common genes found between query and reference data")
    }

    # Pre-compute cell type indices
    query_cell_indices <- split(seq_len(ncol(query_data)),
                                query_data[[query_cell_type_col]])
    ref_cell_indices <- split(seq_len(ncol(reference_data)),
                              reference_data[[ref_cell_type_col]])

    # Filter to only requested cell types
    query_cell_indices <- query_cell_indices[names(query_cell_indices) %in% cell_types]
    ref_cell_indices <- ref_cell_indices[names(ref_cell_indices) %in% cell_types]

    # Get cell types that exist in both datasets
    available_cell_types <- intersect(names(query_cell_indices), names(ref_cell_indices))

    if (length(available_cell_types) == 0) {
        warning("No common cell types found between datasets")
        return(list())
    }

    # Collect all top loading genes across PCs
    all_top_genes <- character(0)
    gene_metadata_list <- list()

    for (pc in pc_subset) {
        pc_loadings <- pca_rotation[common_genes, pc]
        top_loading_indices <- order(abs(pc_loadings),
                                     decreasing = TRUE)[
                                         1:min(n_top_loadings, length(pc_loadings))]
        top_genes <- common_genes[top_loading_indices]
        top_loadings <- pc_loadings[top_loading_indices]

        all_top_genes <- unique(c(all_top_genes, top_genes))

        # Store gene metadata for this PC
        gene_metadata_list[[length(gene_metadata_list) + 1]] <- data.frame(
            gene = top_genes,
            pc = pc,
            loading = top_loadings,
            stringsAsFactors = FALSE
        )
    }

    # Create combined gene metadata
    gene_metadata <- do.call(rbind, gene_metadata_list)

    # Get relevant cells (only from analyzed cell types)
    query_relevant_cells <- unlist(query_cell_indices[available_cell_types],
                                   use.names = FALSE)
    ref_relevant_cells <- unlist(ref_cell_indices[available_cell_types],
                                 use.names = FALSE)

    # Extract expression data for all top genes and relevant cells
    query_expr <- assay(query_data, assay_name)[all_top_genes,
                                                query_relevant_cells, drop = FALSE]
    ref_expr <- assay(reference_data, assay_name)[all_top_genes,
                                                  ref_relevant_cells, drop = FALSE]

    # Create combined expression matrix
    expression_data <- cbind(ref_expr, query_expr)

    # Create cell metadata
    ref_cell_metadata <- data.frame(
        cell_id = colnames(ref_expr),
        dataset = "Reference",
        cell_type = reference_data[[ref_cell_type_col]][ref_relevant_cells],
        original_index = ref_relevant_cells,
        stringsAsFactors = FALSE
    )

    query_cell_metadata <- data.frame(
        cell_id = colnames(query_expr),
        dataset = "Query",
        cell_type = query_data[[query_cell_type_col]][query_relevant_cells],
        original_index = query_relevant_cells,
        stringsAsFactors = FALSE
    )

    cell_metadata <- rbind(ref_cell_metadata, query_cell_metadata)

    # Pre-extract assay data for analysis (using all_top_genes)
    query_assay <- assay(query_data, assay_name)[all_top_genes, , drop = FALSE]
    ref_assay <- assay(reference_data, assay_name)[all_top_genes, , drop = FALSE]

    # Initialize results list for each PC
    pc_results <- vector("list", length(pc_subset))
    names(pc_results) <- paste0("PC", pc_subset)

    # Process each PC
    for (pc_idx in seq_along(pc_subset)) {
        pc <- pc_subset[pc_idx]
        pc_name <- paste0("PC", pc)

        # Get loadings for current PC
        pc_loadings <- pca_rotation[common_genes, pc]

        # Get top loading genes (by absolute value)
        top_loading_indices <- order(abs(pc_loadings),
                                     decreasing = TRUE)[
                                         1:min(n_top_loadings, length(pc_loadings))]
        top_genes <- common_genes[top_loading_indices]
        top_loadings <- pc_loadings[top_loading_indices]

        # Initialize results for this PC
        pc_result_list <- list()

        # Process each cell type
        for (cell_type in available_cell_types) {
            # Get cell indices for current cell type
            query_cells <- query_cell_indices[[cell_type]]
            ref_cells <- ref_cell_indices[[cell_type]]

            # Skip if insufficient cells
            if (length(query_cells) < 3 || length(ref_cells) < 3) {
                next
            }

            # Extract expression data for all top genes at once
            query_expr_matrix <- query_assay[top_genes, query_cells, drop = FALSE]
            ref_expr_matrix <- ref_assay[top_genes, ref_cells, drop = FALSE]

            # Process genes
            gene_results <- processGenesSimple(
                top_genes = top_genes,
                top_loadings = top_loadings,
                query_expr_matrix = query_expr_matrix,
                ref_expr_matrix = ref_expr_matrix,
                cell_type = cell_type
            )

            if (length(gene_results) > 0) {
                pc_result_list <- c(pc_result_list, gene_results)
            }
        }

        # Combine results for this PC
        if (length(pc_result_list) > 0) {
            pc_results_df <- do.call(rbind, pc_result_list)

            # Adjust p-values for multiple testing within this PC
            pc_results_df[["p_adjusted"]] <- p.adjust(pc_results_df[["p_value"]],
                                                      method = adjust_method)

            # Add significance indicator
            pc_results_df[["significant"]] <- pc_results_df[["p_adjusted"]] <= p_value_threshold

            # Sort by adjusted p-value (most significant first)
            pc_results_df <- pc_results_df[order(pc_results_df[["p_adjusted"]]), ]
            rownames(pc_results_df) <- NULL

            pc_results[[pc_name]] <- pc_results_df
        } else {
            pc_results[[pc_name]] <- data.frame()
        }
    }

    # Combine all results into final list
    final_results <- c(
        pc_results,
        list(
            expression_data = expression_data,
            cell_metadata = cell_metadata,
            gene_metadata = gene_metadata,
            percent_var = pc_percent_var
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
#' @details
#' This function processes a set of top loading genes for a specific cell type by:
#' \itemize{
#'   \item Extracting expression values for each gene from query and reference matrices
#'   \item Removing NA values and checking for sufficient sample sizes
#'   \item Calculating mean expression values for both datasets
#'   \item Performing Wilcoxon rank-sum tests to assess distributional differences
#'   \item Returning results in a standardized data frame format
#' }
#'
#' @param top_genes Character vector of gene names to analyze.
#' @param top_loadings Numeric vector of PC loading values corresponding to the genes.
#' @param query_expr_matrix Numeric matrix of expression values for query data (genes × cells).
#' @param ref_expr_matrix Numeric matrix of expression values for reference data (genes × cells).
#' @param cell_type Character string specifying the cell type being analyzed.
#'
#' @keywords internal
#'
#' @return A list of data frames, where each data frame contains results for one gene with columns:
#' \itemize{
#'   \item gene: Gene symbol
#'   \item loading: PC loading value for the gene
#'   \item cell_type: Cell type analyzed
#'   \item p_value: Raw p-value from Wilcoxon rank-sum test
#'   \item mean_query: Mean expression in query data
#'   \item mean_reference: Mean expression in reference data
#' }
#'
#' @author
#' Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Simplified helper function to process genes
processGenesSimple <- function(top_genes, top_loadings,
                               query_expr_matrix, ref_expr_matrix,
                               cell_type) {

    result_list <- list()

    for (i in seq_along(top_genes)) {
        gene <- top_genes[i]
        loading <- top_loadings[i]

        # Extract expression values for current gene
        query_expr <- query_expr_matrix[gene, ]
        ref_expr <- ref_expr_matrix[gene, ]

        # Remove NA values
        query_expr <- query_expr[!is.na(query_expr)]
        ref_expr <- ref_expr[!is.na(ref_expr)]

        # Skip if insufficient data
        if (length(query_expr) < 3 || length(ref_expr) < 3) {
            next
        }

        # Calculate means
        mean_query <- mean(query_expr)
        mean_reference <- mean(ref_expr)

        # Perform Wilcoxon test
        test_result <- suppressWarnings(
            wilcox.test(query_expr, ref_expr, alternative = "two.sided")
        )

        # Store results
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
