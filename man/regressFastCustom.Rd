% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regressPC.R
\name{regressFastCustom}
\alias{regressFastCustom}
\title{Fast Custom Linear Regression for Principal Components}
\usage{
regressFastCustom(pc, indep_var, df)
}
\arguments{
\item{pc}{A character string specifying the principal component column name in the data frame.}

\item{indep_var}{A character string specifying the independent variable specification.
Options include "cell_type", "cell_type * batch", "cell_type * dataset", or other
interaction specifications.}

\item{df}{A data frame containing the principal component scores and categorical predictors.
Must include columns for the specified PC and predictor variables.}
}
\value{
A list containing:
  \item{coefficients}{A data frame with columns \code{coef}, \code{se}, \code{t}, and \code{p.value}
    containing regression coefficients, standard errors, t-statistics, and p-values.}
  \item{r_squared}{A numeric value representing the R-squared of the model.}
}
\description{
Performs efficient linear regression of principal component scores against categorical
predictors (cell types, batches, datasets, or their interactions) using QR decomposition
for numerical stability and computational efficiency.
}
\details{
This function implements a custom linear regression optimized for categorical predictors
commonly used in single-cell RNA sequencing analysis. It uses QR decomposition instead
of normal equations for improved numerical stability and handles rank-deficient design
matrices gracefully. The function supports various model specifications including:
\itemize{
  \item Simple cell type effects: \code{PC ~ cell_type}
  \item Cell type and batch interactions: \code{PC ~ cell_type * batch}
  \item Cell type and dataset interactions: \code{PC ~ cell_type * dataset}
}

The output format is compatible with \code{speedglm} results to maintain consistency
with existing plotting and analysis workflows.
}
\author{
Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
}
\keyword{internal}
