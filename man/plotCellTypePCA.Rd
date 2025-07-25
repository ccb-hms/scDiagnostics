% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCellTypePCA.R
\name{plotCellTypePCA}
\alias{plotCellTypePCA}
\title{Plot Principal Components for Different Cell Types}
\usage{
plotCellTypePCA(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  cell_types = NULL,
  pc_subset = 1:5,
  assay_name = "logcounts",
  lower_facet = c("scatter", "contour", "ellipse", "blank"),
  diagonal_facet = c("ridge", "density", "boxplot"),
  upper_facet = c("blank", "scatter", "contour", "ellipse"),
  max_cells = 2500
)
}
\arguments{
\item{query_data}{A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.}

\item{reference_data}{A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.}

\item{query_cell_type_col}{The column name in the \code{colData} of \code{query_data} that identifies the cell types.}

\item{ref_cell_type_col}{The column name in the \code{colData} of \code{reference_data} that identifies the cell types.}

\item{cell_types}{A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.}

\item{pc_subset}{A numeric vector specifying which principal components to include in the plot. Default is 1:5.}

\item{assay_name}{Name of the assay on which to perform computations. Default is "logcounts".}

\item{lower_facet}{Type of plot to use for the lower panels. Either "scatter" (default), "contour", "ellipse", or "blank".}

\item{diagonal_facet}{Type of plot to use for the diagonal panels. Either "ridge" (default), "density", or "boxplot".}

\item{upper_facet}{Type of plot to use for the upper panels. Either "blank" (default), "scatter", "contour", or "ellipse".}

\item{max_cells}{Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
Default is 2500.}
}
\value{
A ggmatrix object representing a pairs plot of specified principal components for the given cell types and datasets.
}
\description{
This function plots the principal components for different cell types in the query and reference datasets.
}
\details{
This function projects the query dataset onto the principal component space of the reference dataset and then plots the
specified principal components for the specified cell types.
It uses the `projectPCA` function to perform the projection and \code{GGally} to create the pairs plot.
}
