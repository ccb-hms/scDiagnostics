% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGeneExpressionDimred.R
\name{plotGeneExpressionDimred}
\alias{plotGeneExpressionDimred}
\title{Visualize gene expression on a dimensional reduction plot}
\usage{
plotGeneExpressionDimred(
  se_object,
  method = c("TSNE", "UMAP", "PCA"),
  pc_subset = 1:5,
  feature,
  cell_type_col = NULL,
  cell_types = NULL,
  assay_name = "logcounts",
  max_cells = 2500
)
}
\arguments{
\item{se_object}{An object of class \code{\linkS4class{SingleCellExperiment}} containing log-transformed expression matrix and other metadata.
It can be either a reference or query dataset.}

\item{method}{The reduction method to use for visualization. It should be one of the supported methods: "TSNE", "UMAP", or "PCA".}

\item{pc_subset}{An optional vector specifying the principal components (PCs) to include in the plot if method = "PCA".
Default is 1:5.}

\item{feature}{A character string representing the name of the gene or feature to be visualized.}

\item{cell_type_col}{The column name in the \code{colData} of \code{se_object} that identifies the cell types.}

\item{cell_types}{A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.}

\item{assay_name}{Name of the assay on which to perform computations. Default is "logcounts".}

\item{max_cells}{Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
Default is 2500.}
}
\value{
A ggplot object representing the dimensional reduction plot with gene expression.
}
\description{
This function plots gene expression on a dimensional reduction plot using methods like t-SNE, UMAP, or PCA. Each single cell is color-coded based on the expression of a specific gene or feature.
}
\examples{
# Load data
data("query_data")

# Plot gene expression on PCA plot
plotGeneExpressionDimred(se_object = query_data,
                         cell_type_col = "SingleR_annotation",
                         method = "PCA",
                         pc_subset = 1:5,
                         feature = "CD8A",
                         cell_types = "CD4")

}
