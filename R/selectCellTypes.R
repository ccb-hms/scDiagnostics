#' @title Cell Type Selection and Validation for SingleCellExperiment Analysis
#'
#' @description
#' This function selects and validates cell types for functions that analyze \code{\linkS4class{SingleCellExperiment}}
#' objects. It determines which cell types to include based on availability in datasets, applies filtering
#' criteria, and optionally selects the top cell types by cell count.
#'
#' @details
#' The function performs the following selection and validation steps:
#' \itemize{
#'  \item Validates that at least one of \code{query_data} or \code{reference_data} is provided.
#'  \item When \code{dual_only} is TRUE, ensures both datasets are provided.
#'  \item Determines available cell types based on dataset availability and \code{dual_only} setting.
#'  \item If \code{cell_types} is NULL and both datasets are available, includes cell types based on \code{dual_only}.
#'  \item If \code{cell_types} is NULL and only one dataset is available, includes all cell types from that dataset.
#'  \item If \code{cell_types} is provided, filters to include only valid types.
#'  \item If \code{n_cell_types} is specified, selects the top cell types by total cell count.
#'  \item Returns the selected and validated cell types as character strings.
#' }
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' Can be \code{NULL} if only reference data is available.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' Can be \code{NULL} if only query data is available.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types. Should be \code{NULL} if \code{query_data} is \code{NULL}.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types. Should be \code{NULL} if \code{reference_data} is \code{NULL}.
#' @param cell_types A character vector specifying the cell types to validate. If \code{NULL},
#' cell types will be automatically selected based on dataset availability and \code{dual_only} setting.
#' @param dual_only A logical value indicating whether cell types must be present in both datasets.
#' If \code{TRUE}, both \code{query_data} and \code{reference_data} must be provided, and only
#' cell types present in both datasets will be considered. Default is \code{FALSE}.
#' @param n_cell_types An integer specifying the maximum number of cell types to select based on
#' highest cell count. If \code{NULL}, all valid cell types are returned. Default is \code{NULL}.
#'
#' @keywords internal
#'
#' @return A character vector of selected cell types that meet the specified criteria.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to select and validate cell types
selectCellTypes <- function(query_data = NULL,
                            reference_data = NULL,
                            query_cell_type_col = NULL,
                            ref_cell_type_col = NULL,
                            cell_types = NULL,
                            dual_only = FALSE,
                            n_cell_types = NULL) {

    # Check that at least one dataset is provided
    if (is.null(query_data) && is.null(reference_data)) {
        stop("At least one of 'query_data' or 'reference_data' must be provided.")
    }

    # Check dual_only requirements
    if (dual_only && (is.null(query_data) || is.null(reference_data))) {
        stop("When 'dual_only' is TRUE, both 'query_data' and 'reference_data' must be provided.")
    }

    # Validate n_cell_types if provided
    if (!is.null(n_cell_types)) {
        if (!is.numeric(n_cell_types) || length(n_cell_types) != 1 ||
            n_cell_types <= 0 || n_cell_types != as.integer(n_cell_types)) {
            stop("'n_cell_types' must be a positive integer.")
        }
    }

    # Extract available cell types from each dataset and convert to character
    available_query_types <- NULL
    available_ref_types <- NULL

    if (!is.null(query_data) && !is.null(query_cell_type_col)) {
        if (!query_cell_type_col %in% names(colData(query_data))) {
            stop("'query_cell_type_col' is not found in query_data colData.")
        }
        available_query_types <- as.character(na.omit(unique(query_data[[query_cell_type_col]])))
    }

    if (!is.null(reference_data) && !is.null(ref_cell_type_col)) {
        if (!ref_cell_type_col %in% names(colData(reference_data))) {
            stop("'ref_cell_type_col' is not found in reference_data colData.")
        }
        available_ref_types <- as.character(na.omit(unique(reference_data[[ref_cell_type_col]])))
    }

    # Determine valid cell types based on input and availability
    if (is.null(cell_types)) {
        # Auto-select cell types based on dataset availability
        if (!is.null(available_query_types) && !is.null(available_ref_types)) {
            # Both datasets available - use dual_only setting
            if (dual_only) {
                all_available_types <- intersect(available_query_types, available_ref_types)
            } else {
                all_available_types <- unique(c(available_query_types, available_ref_types))
            }
        } else {
            # Only one dataset available - use all cell types from that dataset
            all_available_types <- unique(c(available_query_types, available_ref_types))
        }
    } else {
        # User provided cell types - convert to character and validate against availability
        cell_types <- as.character(cell_types)

        # Determine which types are valid based on dual_only setting
        if (dual_only) {
            if (is.null(available_query_types) || is.null(available_ref_types)) {
                stop("When 'dual_only' is TRUE, both datasets must have valid cell type columns.")
            }
            valid_types_pool <- intersect(available_query_types, available_ref_types)
        } else {
            valid_types_pool <- unique(c(available_query_types, available_ref_types))
        }

        # Filter to only those present in valid pool
        all_available_types <- cell_types[cell_types %in% valid_types_pool]

        # Check if any requested cell types were found
        if (length(all_available_types) == 0) {
            missing_types <- setdiff(cell_types, valid_types_pool)
            if (dual_only) {
                stop("None of the specified cell types are present in both datasets. ",
                     "Missing types: ", paste(missing_types, collapse = ", "))
            } else {
                stop("None of the specified cell types are present in the provided datasets. ",
                     "Missing types: ", paste(missing_types, collapse = ", "))
            }
        }

        # Warn about missing cell types if some but not all were found
        missing_types <- setdiff(cell_types, all_available_types)
        if (length(missing_types) > 0) {
            if (dual_only) {
                warning("Some specified cell types are not present in both datasets and will be excluded: ",
                        paste(missing_types, collapse = ", "))
            } else {
                warning("Some specified cell types are not present in the provided datasets and will be excluded: ",
                        paste(missing_types, collapse = ", "))
            }
        }
    }

    # Ensure all_available_types is character and remove any NA values
    all_available_types <- as.character(all_available_types)
    all_available_types <- all_available_types[!is.na(all_available_types)]

    # If no cell types are available, return error
    if (length(all_available_types) == 0) {
        if (dual_only) {
            stop("No common cell types found between query and reference datasets.")
        } else {
            stop("No valid cell types found in the provided datasets.")
        }
    }

    # Apply n_cell_types selection if specified
    if (!is.null(n_cell_types) && length(all_available_types) > n_cell_types) {
        # Count cells for each cell type
        cell_counts <- numeric(length(all_available_types))
        names(cell_counts) <- all_available_types

        for (ct in all_available_types) {
            count <- 0

            # Count cells in query data if available and cell type is present
            if (!is.null(available_query_types) && ct %in% available_query_types) {
                count <- count + sum(as.character(query_data[[query_cell_type_col]]) == ct, na.rm = TRUE)
            }

            # Count cells in reference data if available and cell type is present
            if (!is.null(available_ref_types) && ct %in% available_ref_types) {
                count <- count + sum(as.character(reference_data[[ref_cell_type_col]]) == ct, na.rm = TRUE)
            }

            cell_counts[ct] <- count
        }

        # Select top n_cell_types by count
        top_indices <- order(cell_counts, decreasing = TRUE)[1:n_cell_types]
        selected_types <- names(cell_counts)[top_indices]

        # Warn about selection
        excluded_types <- setdiff(all_available_types, selected_types)
        warning("Selected top ", n_cell_types, " cell types by cell count. ",
                "Excluded types: ", paste(excluded_types, collapse = ", "))

        final_cell_types <- selected_types
    } else {
        final_cell_types <- all_available_types
    }

    # Ensure final result is character and sort for consistent output
    final_cell_types <- as.character(final_cell_types)
    final_cell_types <- sort(final_cell_types)

    return(final_cell_types)
}
