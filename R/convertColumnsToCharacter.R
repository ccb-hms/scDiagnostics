#' @title Convert Specified Columns to Character in SingleCellExperiment Objects
#'
#' @description
#' This function converts specified columns in the \code{colData} of a \code{\linkS4class{SingleCellExperiment}}
#' object to character type. It checks that the specified columns exist and only performs conversion
#' when necessary (i.e., when columns are not already character type).
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'  \item Validates that the input is a \code{\linkS4class{SingleCellExperiment}} object.
#'  \item Checks that all specified columns exist in the \code{colData} of the object.
#'  \item Converts each specified column to character type if it is not already character.
#'  \item Returns the modified \code{\linkS4class{SingleCellExperiment}} object.
#'  \item If all specified columns are already character type, returns the object unchanged.
#' }
#'
#' This function is particularly useful for handling factor columns that need to be converted
#' to character for downstream analysis functions that expect character input.
#'
#' @param sce_object A \code{\linkS4class{SingleCellExperiment}} object containing single-cell data.
#' @param convert_cols A character vector specifying the column names in \code{colData} to convert
#' to character type. All specified columns must exist in the \code{colData}.
#'
#' @keywords internal
#'
#' @return A \code{\linkS4class{SingleCellExperiment}} object with the specified columns converted
#' to character type in the \code{colData}.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to convert specified columns to character in SingleCellExperiment objects
convertColumnsToCharacter <- function(sce_object, convert_cols) {

    # Validate input object
    if (!is(sce_object, "SingleCellExperiment")) {
        stop("'sce_object' must be a SingleCellExperiment object.")
    }

    # Validate convert_cols parameter
    if (!is.character(convert_cols) || length(convert_cols) == 0) {
        stop("'convert_cols' must be a non-empty character vector.")
    }

    # Check if all specified columns exist in colData
    existing_cols <- names(SummarizedExperiment::colData(sce_object))
    missing_cols <- convert_cols[!convert_cols %in% existing_cols]

    if (length(missing_cols) > 0) {
        stop("The following columns are not found in colData: ",
             paste(missing_cols, collapse = ", "))
    }

    # Convert columns to character if they are not already character
    col_data <- SummarizedExperiment::colData(sce_object)
    conversion_needed <- FALSE

    for (col in convert_cols) {
        if (!is.character(col_data[[col]])) {
            col_data[[col]] <- as.character(col_data[[col]])
            conversion_needed <- TRUE
        }
    }

    # Update the object only if conversion was needed
    if (conversion_needed) {
        SummarizedExperiment::colData(sce_object) <- col_data
    }

    return(sce_object)
}
