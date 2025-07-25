% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateAveragePairwiseCorrelation.R,
%   R/plot.calculateAveragePairwiseCorrelationObject.R
\name{calculateAveragePairwiseCorrelation}
\alias{calculateAveragePairwiseCorrelation}
\alias{plot.calculateAveragePairwiseCorrelationObject}
\title{Compute Average Pairwise Correlation between Cell Types}
\usage{
calculateAveragePairwiseCorrelation(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  pc_subset = 1:10,
  correlation_method = c("spearman", "pearson"),
  assay_name = "logcounts",
  max_cells = 2500
)

\method{plot}{calculateAveragePairwiseCorrelationObject}(x, ...)
}
\arguments{
\item{query_data}{A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.}

\item{reference_data}{A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.}

\item{query_cell_type_col}{The column name in the \code{colData} of \code{query_data} that identifies the cell types.}

\item{ref_cell_type_col}{The column name in the \code{colData} of \code{reference_data} that identifies the cell types.}

\item{cell_types}{A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.}

\item{pc_subset}{A numeric vector specifying which principal components to use in the analysis. Default is 1:10.
If set to \code{NULL} then no dimensionality reduction is performed and the assay data is used directly for computations.}

\item{correlation_method}{The correlation method to use for calculating pairwise correlations.}

\item{assay_name}{Name of the assay on which to perform computations. Default is "logcounts".}

\item{max_cells}{Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
Default is 2500.}

\item{x}{Output matrix from `calculateAveragePairwiseCorrelation` function.}

\item{...}{Additional arguments to be passed to the plotting function.}
}
\value{
A matrix containing the average pairwise correlation values.
        Rows and columns are labeled with the cell types. Each element
        in the matrix represents the average correlation between a pair
        of cell types.

The S3 plot method returns a \code{ggplot} object representing the heatmap plot.
}
\description{
Computes the average pairwise correlations between specified cell types
in single-cell gene expression data.

The S3 plot method takes the output of the `calculateAveragePairwiseCorrelation` function,
which should be a matrix of pairwise correlations, and plots it as a heatmap.
}
\details{
This function operates on \code{\linkS4class{SingleCellExperiment}} objects,
ideal for single-cell analysis workflows. It calculates pairwise correlations between query and
reference cells using a specified correlation method, then averages these correlations for each
cell type pair. This function aids in assessing the similarity between cells in reference and query datasets,
providing insights into the reliability of cell type annotations in single-cell gene expression data.

The S3 plot method converts the correlation matrix into a dataframe, creates a heatmap using ggplot2,
and customizes the appearance of the heatmap with updated colors and improved aesthetics.
}
\examples{
# Load data
data("reference_data")
data("query_data")

# Compute pairwise correlations
cor_matrix_avg <- calculateAveragePairwiseCorrelation(query_data = query_data,
                                                      reference_data = reference_data,
                                                      query_cell_type_col = "SingleR_annotation",
                                                      ref_cell_type_col = "expert_annotation",
                                                      cell_types = c("CD4", "CD8", "B_and_plasma"),
                                                      pc_subset = 1:10,
                                                      correlation_method = "spearman")

# Visualize correlation output
plot(cor_matrix_avg)

}
\seealso{
\code{\link{plot.calculateAveragePairwiseCorrelationObject}}

\code{\link{calculateAveragePairwiseCorrelation}}
}
\author{
Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
}
\keyword{internal}
