#' @title Reference Single-Cell RNA-Seq Dataset
#'
#' @description
#' This dataset contains the processed reference dataset from the HeOrganAtlas dataset for Marrow tissue.
#' It has been preprocessed to include log-normalized counts, specific metadata columns, and PCA, t-SNE, and UMAP results.
#'
#' @details
#' This dataset underwent the following steps:
#' \itemize{
#'   \item Loads the HeOrganAtlas dataset specifically for Marrow tissue from the \code{scRNAseq} package.
#'   \item Divides the loaded dataset into a reference dataset used for downstream analysis.
#'   \item Performs log normalization on the reference dataset using the function \code{logNormCounts} from the \code{scuttle} package.
#'   \item Selects the column \code{expert_annotation}) from the cell metadata for downstream analysis.
#'   \item Selects highly variable genes (HVGs) using the function \code{getTopHVGs} from the \code{scran} package on the reference dataset.
#'   \item Performs Principal Component Analysis (PCA) on the reference dataset using the function \code{runPCA} from the \code{scater} package.
#'   \item Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on the reference dataset using the function \code{runTSNE} from the \code{scater} package.
#'   \item Performs Uniform Manifold Approximation and Projection (UMAP) on the reference dataset using the function \code{runUMAP} from the \code{scater} package.
#' }
#'
#' @seealso Use \code{data("reference_data")} to load and access the resulting reference dataset.
#'
#' @source The HeOrganAtlas dataset, available through the scRNAseq package.

#' @references He, et al. (2020). HeOrganAtlas: a comprehensive human organ atlas based on single-cell RNA sequencing.
#'
#' @examples
#' # Load and explore the reference dataset
#' data("reference_data")
#'
#' @keywords internal
#'
"reference_data"

#' @title Query Single-Cell RNA-Seq Dataset
#'
#' @description
#' This dataset contains the processed query dataset from the HeOrganAtlas dataset for Marrow tissue.
#' It has been preprocessed to include log-normalized counts, specific metadata columns, annotations
#' based on SingleR cell type scoring, and PCA, t-SNE, and UMAP results.
#'
#' @details
#' This dataset underwent the following steps:
#' \itemize{
#'   \item Loads the HeOrganAtlas dataset specifically for Marrow tissue from the \code{scRNAseq} package.
#'   \item Divides the loaded dataset into a query dataset used for downstream analysis.
#'   \item Performs log normalization on the query dataset using the function \code{logNormCounts} from the \code{scuttle} package.
#'   \item Selects specific columns (\code{percent_mito}, \code{expert_annotation}) from the cell metadata for downstream analysis.
#'   \item Selects highly variable genes (HVGs) using the function \code{getTopHVGs} from the \code{scran} package on the query dataset.
#'   \item Computes AUC gene set scores using the function \code{AUCell_calcAUC} from the \code{AUCell} package based on a CD4 T cell
#'         signature containing 12 known CD4 T cell marker genes (IL7R, CCR7, SELL, LEF1, TCF7, LTB, KLF2, IL32, CD2, CD3D, CD3E, CD3G)
#'         and adds these scores to the query dataset as \code{gene_set_scores}.
#'   \item Intersects the highly variable genes between the query and reference datasets to obtain common genes for analysis.
#'   \item Performs Principal Component Analysis (PCA) on the query dataset using the function \code{runPCA} from the \code{scater}
#'         package.
#'   \item Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on the query dataset using the function \code{runTSNE} from the
#'         \code{scater} package.
#'   \item Performs Uniform Manifold Approximation and Projection (UMAP) on the query dataset using the function \code{runUMAP} from the
#'         \code{scater} package.
#'   \item Adds SingleR annotations (\code{SingleR_annotation}) and annotation scores (\code{annotation_scores}) to the query dataset
#'         using the function \code{SingleR} from the \code{SingleR} package.

#' }
#'
#' @seealso
#' Use \code{data("query_data")} to load and access the resulting query dataset and the
#' \code{data("reference_data")} for comparison with the reference dataset.
#'
#' @source The HeOrganAtlas dataset, available through the scRNAseq package.
#' @references He, et al. (2020). HeOrganAtlas: a comprehensive human organ atlas based on single-cell RNA sequencing.
#'
#' @examples
#' # Load and explore the query dataset
#' data("query_data")
#'
#' @keywords internal
#'
"query_data"

#' @title Quality Control Single-Cell RNA-Seq Dataset
#'
#' @description
#' This dataset contains the processed query dataset from the Bunis haematopoietic stem and progenitor cell data.
#' It has been preprocessed to include log-normalized counts, QC metrics, SingleR cell type predictions,
#' and annotation scores.
#'
#' @details
#' This dataset underwent the following steps:
#' \itemize{
#'   \item Loads the \code{hpca} reference dataset using \code{fetchReference} from the \code{celldex} package.
#'   \item Loads the QC dataset (Bunis haematopoietic stem and progenitor cell data) from Bunis DG et al. (2021).
#'   \item Adds QC metrics to the QC dataset using the function \code{addPerCellQCMetrics} from the \code{scuttle} package.
#'   \item Performs log normalization on the QC dataset using the function \code{logNormCounts} from the \code{scuttle} package.
#'   \item Runs SingleR to predict cell types and assigns predicted labels to the QC dataset using the function \code{SingleR} from the \code{SingleR} package.
#'   \item Assigns annotation scores to the QC dataset.
#'   \item Selects specific columns (\code{total}, \code{SingleR_annotation}, \code{annotation_scores}) from the cell metadata for downstream analysis.
#'   \item Selects highly variable genes (HVGs) using the function \code{getTopHVGs} from the \code{scran} package on the QC dataset.
#' }
#'
#' @seealso
#' Use \code{data("qc_data")} to load and access the resulting quality control dataset.
#'
#' @source Bunis DG et al. (2021). Single-Cell Mapping of Progressive Fetal-to-Adult Transition in Human Naive T Cells Cell Rep. 34(1): 108573
#' @references Bunis DG et al. (2021). Single-Cell Mapping of Progressive Fetal-to-Adult Transition in Human Naive T Cells Cell Rep. 34(1): 108573
#'
#' @examples
#' # Load and explore the quality control dataset
#' data("qc_data")
#'
#' @keywords internal
#'
"qc_data"
