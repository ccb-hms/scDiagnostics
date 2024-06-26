% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateSampleDistancesSimilarity.R
\name{calculateSampleDistancesSimilarity}
\alias{calculateSampleDistancesSimilarity}
\title{Function to Calculate Bhattacharyya Coefficients and Hellinger Distances}
\usage{
calculateSampleDistancesSimilarity(
  query_data,
  reference_data,
  query_cell_type_col,
  ref_cell_type_col,
  sample_names,
  n_components = 10,
  pc_subset = c(1:5)
)
}
\arguments{
\item{query_data}{A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.}

\item{reference_data}{A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.}

\item{query_cell_type_col}{character. The column name in the \code{colData} of \code{query_data} 
that identifies the cell types.}

\item{ref_cell_type_col}{character. The column name in the \code{colData} of \code{reference_data} 
that identifies the cell types.}

\item{sample_names}{A character vector specifying the names of the query samples for which to compute distance measures.}

\item{n_components}{An integer specifying the number of principal components to use for projection. Defaults to 10.}

\item{pc_subset}{A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.}
}
\value{
A list containing distance data for each cell type. Each entry in the list contains:
\describe{
  \item{ref_distances}{A vector of all pairwise distances within the reference subset for the cell type.}
  \item{query_to_ref_distances}{A matrix of distances from each query sample to all reference samples for the cell type.}
}
}
\description{
This function computes Bhattacharyya coefficients and Hellinger distances to quantify the similarity of density 
distributions between query samples and reference data for each cell type.
}
\details{
This function first computes distance data using the \code{calculateSampleDistances} function, which calculates 
pairwise distances between samples within the reference data and between query samples and reference samples in the PCA space.
Bhattacharyya coefficients and Hellinger distances are calculated to quantify the similarity of density distributions between query 
samples and reference data for each cell type. Bhattacharyya coefficient measures the similarity of two probability distributions, 
while Hellinger distance measures the distance between two probability distributions.

Bhattacharyya coefficients range between 0 and 1. A value closer to 1 indicates higher similarity between distributions, while a value 
closer to 0 indicates lower similarity

Hellinger distances range between 0 and 1. A value closer to 0 indicates higher similarity between distributions, while a value 
closer to 1 indicates lower similarity.
}
\examples{
# Load required libraries
library(scRNAseq)
library(scuttle)
library(SingleR)
library(scran)
library(scater)

# Load data (replace with your data loading)
sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)

# Divide the data into reference and query datasets
set.seed(100)
indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
ref_data <- sce[, indices]
query_data <- sce[, -indices]

# log transform datasets
ref_data <- scuttle::logNormCounts(ref_data)
query_data <- scuttle::logNormCounts(query_data)

# Get cell type scores using SingleR (or any other cell type annotation method)
scores <- SingleR::SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)

# Add labels to query object
colData(query_data)$labels <- scores$labels

# Selecting highly variable genes (can be customized by the user)
ref_var <- scran::getTopHVGs(ref_data, n = 2000)
query_var <- scran::getTopHVGs(query_data, n = 2000)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_var, query_var)
ref_data_subset <- ref_data[common_genes, ]
query_data_subset <- query_data[common_genes, ]

# Run PCA on the reference data
ref_data_subset <- runPCA(ref_data_subset, ncomponents = 50)

# Plot the PC data
distance_data <- calculateSampleDistances(query_data_subset, ref_data_subset, 
                                          n_components = 10, 
                                          query_cell_type_col = "labels", 
                                          ref_cell_type_col = "reclustered.broad",
                                          pc_subset = c(1:10)) 

# Identify outliers for CD4
cd4_anomalies <- detectAnomaly(ref_data_subset, query_data_subset, 
                               query_cell_type_col = "labels", 
                               ref_cell_type_col = "reclustered.broad",
                               n_components = 10,
                               n_tree = 500,
                               anomaly_treshold = 0.5)$CD4
cd4_top5_anomalies <- names(sort(cd4_anomalies$query_anomaly_scores, decreasing = TRUE)[1:6])

# Get overlap measures
overlap_measures <- calculateSampleDistancesSimilarity(query_data_subset,ref_data_subset, 
                                                       sample_names = cd4_top5_anomalies,
                                                       n_components = 10, 
                                                       query_cell_type_col = "labels", 
                                                       ref_cell_type_col = "reclustered.broad",
                                                       pc_subset = c(1:10))
overlap_measures


}
\author{
Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
}
