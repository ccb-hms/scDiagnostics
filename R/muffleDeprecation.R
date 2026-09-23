#' @title Evaluate an Expression Without Dependency Deprecation Warnings
#'
#' @description This function evaluates an expression with deprecation
#' warnings muffled, leaving every other condition untouched.
#'
#' @details The highly variable gene helpers of the \code{scran} package
#' (\code{modelGeneVar()}, \code{getTopHVGs()}, and the \code{fitTrendVar()}
#' that \code{modelGeneVar()} calls internally) are deprecated in favour of
#' the \code{scrapper} package as of \code{scran} 1.41. They still work, but
#' every call signals a condition of class \code{deprecatedWarning}, which
#' \code{R CMD check} reports as a significant warning in the examples. Only
#' that condition class is muffled here, so genuine warnings still reach the
#' user. This helper becomes unnecessary once the highly variable gene
#' selection is ported to \code{scrapper}.
#'
#' @param expr An expression to evaluate.
#'
#' @keywords internal
#'
#' @return The value of \code{expr}.
#'
#' @author Anthony Christidis,
#' \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @noRd
#'
# Function to evaluate an expression without deprecation warnings
muffleDeprecation <- function(expr) {
    withCallingHandlers(
        expr,
        deprecatedWarning = function(w) invokeRestart("muffleWarning")
    )
}
