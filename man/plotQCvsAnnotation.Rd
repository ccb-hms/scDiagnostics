% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotQCvsAnnotation.R
\name{plotQCvsAnnotation}
\alias{plotQCvsAnnotation}
\title{Scatter plot: QC stats vs Cell Type Annotation Scores}
\usage{
plotQCvsAnnotation(
  se_object,
  cell_type_col,
  cell_types = NULL,
  qc_col,
  score_col,
  max_cells = 2500
)
}
\arguments{
\item{se_object}{A \code{\linkS4class{SingleCellExperiment}} containing the single-cell
expression data and metadata.}

\item{cell_type_col}{The column name in the \code{colData} of \code{se_object}
that contains the cell type labels.}

\item{cell_types}{A vector of cell type labels to plot (e.g., c("T-cell", "B-cell")).
Defaults to \code{NULL}, which will include all the cells.}

\item{qc_col}{A column name in the \code{colData} of \code{se_object} that
contains the QC stats of interest.}

\item{score_col}{The column name in the \code{colData} of \code{se_object} that
contains the cell type annotation scores.}

\item{max_cells}{Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
Default is 2500.}
}
\value{
A ggplot object displaying a scatter plot of QC stats vs annotation scores,
        where each point represents a cell, color-coded by its cell type.
}
\description{
Creates a scatter plot to visualize the relationship between QC stats (e.g., library size)
and cell type annotation scores for one or more cell types.
}
\details{
This function generates a scatter plot to explore the relationship between various quality
control (QC) statistics, such as library size and mitochondrial percentage, and cell type
annotation scores. By examining these relationships, users can assess whether specific QC
metrics, systematically influence the confidence in cell type annotations,
which is essential for ensuring reliable cell type annotation.
}
\examples{
# Load data
data("qc_data")

# Remove cell types with very few cells
qc_data_subset <- qc_data[, !(qc_data$SingleR_annotation
                              \%in\% c("Chondrocytes", "DC",
                                     "Neurons","Platelets"))]

p1 <- plotQCvsAnnotation(se_object = qc_data_subset,
                         cell_type_col = "SingleR_annotation",
                         cell_types = NULL,
                         qc_col = "total",
                         score_col = "annotation_scores")
p1 + ggplot2::xlab("Library Size")

}
