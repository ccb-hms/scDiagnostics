#' @title Reference Single-Cell RNA-Seq Dataset
#'
#' @description This dataset contains the processed reference dataset from the
#' HeOrganAtlas dataset for Marrow tissue. It has been preprocessed to include
#' log-normalized counts, specific metadata columns, and PCA, t-SNE, and UMAP
#' results.
#'
#' @details This dataset underwent the following steps:
#' \itemize{
#'  \item Loads the HeOrganAtlas dataset specifically for Marrow tissue from the
#'  \code{scRNAseq} package.
#'  \item Divides the loaded dataset into a reference dataset used for
#'  downstream analysis.
#'  \item Performs log normalization on the reference dataset using the function
#'  \code{logNormCounts} from the \code{scuttle} package.
#'  \item Selects the column \code{expert_annotation}) from the cell metadata
#'  for downstream analysis.
#'  \item Selects highly variable genes (HVGs) using the function
#'  \code{getTopHVGs} from the \code{scran} package on the reference dataset.
#'  \item Performs Principal Component Analysis (PCA) on the reference dataset
#'  using the function \code{runPCA} from the \code{scater} package.
#'  \item Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on the
#'  reference dataset using the function \code{runTSNE} from the \code{scater}
#'  package.
#'  \item Performs Uniform Manifold Approximation and Projection (UMAP) on the
#'  reference dataset using the function \code{runUMAP} from the \code{scater}
#'  package.
#' }
#'
#' @seealso Use \code{data("reference_data")} to load and access the resulting
#' reference dataset.
#'
#' @source The HeOrganAtlas dataset, available through the scRNAseq package.

#' @references He, et al. (2020). HeOrganAtlas: a comprehensive human organ
#' atlas based on single-cell RNA sequencing.
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
#' @description This dataset contains the processed query dataset from the
#' HeOrganAtlas dataset for Marrow tissue. It has been preprocessed to include
#' log-normalized counts, specific metadata columns, annotations based on
#' SingleR cell type scoring, and PCA, t-SNE, and UMAP results.
#'
#' @details This dataset underwent the following steps:
#' \itemize{
#'  \item Loads the HeOrganAtlas dataset specifically for Marrow tissue from the
#'  \code{scRNAseq} package.
#'  \item Divides the loaded dataset into a query dataset used for downstream
#'  analysis.
#'  \item Performs log normalization on the query dataset using the function
#'  \code{logNormCounts} from the \code{scuttle} package.
#'  \item Selects specific columns (\code{percent_mito},
#'  \code{expert_annotation}) from the cell metadata for downstream analysis.
#'  \item Selects highly variable genes (HVGs) using the function
#'  \code{getTopHVGs} from the \code{scran} package on the query dataset.
#'  \item Computes AUC gene set scores using the function \code{AUCell_calcAUC}
#'  from the \code{AUCell} package based on a CD4 T cell signature containing 12
#'  known CD4 T cell marker genes (IL7R, CCR7, SELL, LEF1, TCF7, LTB, KLF2,
#'  IL32, CD2, CD3D, CD3E, CD3G) and adds these scores to the query dataset as
#'  \code{gene_set_scores}.
#'  \item Intersects the highly variable genes between the query and reference
#'  datasets to obtain common genes for analysis.
#'  \item Performs Principal Component Analysis (PCA) on the query dataset using
#'  the function \code{runPCA} from the \code{scater} package.
#'  \item Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on the
#'  query dataset using the function \code{runTSNE} from the \code{scater}
#'  package.
#'  \item Performs Uniform Manifold Approximation and Projection (UMAP) on the
#'  query dataset using the function \code{runUMAP} from the \code{scater}
#'  package.
#'  \item Adds SingleR annotations (\code{SingleR_annotation}) and annotation
#'  scores (\code{annotation_scores}) to the query dataset using the function
#'  \code{SingleR} from the \code{SingleR} package.

#' }
#'
#' @seealso Use \code{data("query_data")} to load and access the resulting query
#' dataset and the \code{data("reference_data")} for comparison with the
#' reference dataset.
#'
#' @source The HeOrganAtlas dataset, available through the scRNAseq package.
#' @references He, et al. (2020). HeOrganAtlas: a comprehensive human organ
#' atlas based on single-cell RNA sequencing.
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
#' @description This dataset contains the processed query dataset from the Bunis
#' haematopoietic stem and progenitor cell data. It has been preprocessed to
#' include log-normalized counts, QC metrics, SingleR cell type predictions, and
#' annotation scores.
#'
#' @details This dataset underwent the following steps:
#' \itemize{
#'  \item Loads the \code{hpca} reference dataset using \code{fetchReference}
#'  from the \code{celldex} package.
#'  \item Loads the QC dataset (Bunis haematopoietic stem and progenitor cell
#'  data) from Bunis DG et al. (2021).
#'  \item Adds QC metrics to the QC dataset using the function
#'  \code{addPerCellQCMetrics} from the \code{scuttle} package.
#'  \item Performs log normalization on the QC dataset using the function
#'  \code{logNormCounts} from the \code{scuttle} package.
#'  \item Runs SingleR to predict cell types and assigns predicted labels to the
#'  QC dataset using the function \code{SingleR} from the \code{SingleR}
#'  package.
#'  \item Assigns annotation scores to the QC dataset.
#'  \item Selects specific columns (\code{total}, \code{SingleR_annotation},
#'  \code{annotation_scores}) from the cell metadata for downstream analysis.
#'  \item Selects highly variable genes (HVGs) using the function
#'  \code{getTopHVGs} from the \code{scran} package on the QC dataset.
#' }
#'
#' @seealso Use \code{data("qc_data")} to load and access the resulting quality
#' control dataset.
#'
#' @source Bunis DG et al. (2021). Single-Cell Mapping of Progressive
#' Fetal-to-Adult Transition in Human Naive T Cells Cell Rep. 34(1): 108573
#' @references Bunis DG et al. (2021). Single-Cell Mapping of Progressive
#' Fetal-to-Adult Transition in Human Naive T Cells Cell Rep. 34(1): 108573
#'
#' @examples
#' # Load and explore the quality control dataset
#' data("qc_data")
#'
#' @keywords internal
#'
"qc_data"

#' @title Zeisel Brain Reference Single-Cell RNA-Seq Dataset
#'
#' @description This dataset contains the processed reference dataset from
#' the Zeisel mouse brain data. It has been preprocessed to include
#' log-normalized counts, a cell type column, and PCA results, and is used to
#' benchmark \code{detectAnomaly} and \code{calculateReconstructionError}
#' against known ground truth.
#'
#' @details This dataset underwent the following steps:
#' \itemize{
#'  \item Loads the Zeisel mouse brain dataset from the \code{scRNAseq}
#'  package.
#'  \item Performs log normalization using the function \code{logNormCounts}
#'  from the \code{scuttle} package.
#'  \item Renames the \code{level1class} column to \code{true_cell_type}.
#'  \item Divides the data into a 70% reference and 30% query split.
#'  \item Selects highly variable genes (HVGs) using the function
#'  \code{getTopHVGs} from the \code{scran} package, intersected between the
#'  reference and query datasets.
#'  \item Performs Principal Component Analysis (PCA) on the reference
#'  dataset using the function \code{runPCA} from the \code{scater} package.
#' }
#'
#' @seealso Use \code{data("zeisel_reference_data")} to load and access the
#' resulting reference dataset, and \code{data("zeisel_query_data")} for the
#' corresponding query dataset.
#'
#' @source The Zeisel mouse brain dataset, available through the scRNAseq
#' package.
#' @references Zeisel A, et al. (2015). Cell types in the mouse cortex and
#' hippocampus revealed by single-cell RNA-seq. Science 347(6226):1138-42.
#'
#' @examples
#' # Load and explore the Zeisel reference dataset
#' data("zeisel_reference_data")
#'
#' @keywords internal
#'
"zeisel_reference_data"

#' @title Zeisel Brain Query Single-Cell RNA-Seq Dataset
#'
#' @description This dataset contains the processed query dataset from the
#' Zeisel mouse brain data, for use alongside \code{zeisel_reference_data}.
#'
#' @details See \code{zeisel_reference_data} for the processing steps shared
#' by both datasets.
#'
#' @seealso Use \code{data("zeisel_query_data")} to load and access the
#' resulting query dataset, and \code{data("zeisel_reference_data")} for
#' comparison with the reference dataset.
#'
#' @source The Zeisel mouse brain dataset, available through the scRNAseq
#' package.
#' @references Zeisel A, et al. (2015). Cell types in the mouse cortex and
#' hippocampus revealed by single-cell RNA-seq. Science 347(6226):1138-42.
#'
#' @examples
#' # Load and explore the Zeisel query dataset
#' data("zeisel_query_data")
#'
#' @keywords internal
#'
"zeisel_query_data"

#' @title Zeisel Brain Anomaly Detection Benchmark Results
#'
#' @description A precomputed benchmark of \code{detectAnomaly} (Isolation
#' Forest) and \code{calculateReconstructionError} on the full Zeisel mouse
#' brain dataset, evaluating detection accuracy against known ground truth
#' (a withheld cell type) across label-noise, class-imbalance, and
#' batch-effect gradients, plus a hyperparameter grid search. Used in the
#' \code{ZeiselBenchmarking} vignette to illustrate how these functions
#' perform beyond a single worked example.
#'
#' @details A named list with three data frames:
#' \itemize{
#'  \item \code{gradients}: AUROC, AUPRC, sensitivity, and specificity for
#'  both methods across 3 baseline scenarios (a distinct, a related, and a
#'  rare cell type withheld from the reference) and, for the related and
#'  rare scenarios, across gradients of label noise (0-40% shuffled labels),
#'  reference class imbalance (300 to 10 cells), and query batch effects
#'  (a mean shift applied to 20% of genes).
#'  \item \code{if_tuning}: sensitivity/specificity of \code{detectAnomaly}
#'  across a grid of PC subsets, HVG counts, and anomaly thresholds.
#'  \item \code{re_tuning}: sensitivity/specificity of
#'  \code{calculateReconstructionError} across a grid of HVG counts, PC
#'  subsets, and MAD thresholds.
#' }
#' Computed once, offline, on the full (non-downsampled) Zeisel dataset; see
#' \code{inst/script/ZeiselBenchmarkResults.R} for the exact procedure.
#'
#' @seealso Use \code{data("zeisel_benchmark_results")} to load and access
#' the benchmark results.
#'
#' @source Computed from the Zeisel mouse brain dataset (scRNAseq package)
#' using \code{scDiagnostics::detectAnomaly} and
#' \code{scDiagnostics::calculateReconstructionError}.
#' @references Zeisel A, et al. (2015). Cell types in the mouse cortex and
#' hippocampus revealed by single-cell RNA-seq. Science 347(6226):1138-42.
#'
#' @examples
#' # Load and explore the Zeisel benchmark results
#' data("zeisel_benchmark_results")
#'
#' @keywords internal
#'
"zeisel_benchmark_results"

#' @title COVID-19 PBMC Reference (Healthy) Single-Cell RNA-Seq Dataset
#'
#' @description This dataset contains a downsampled reference (healthy
#' donor) subset of the Stephenson et al. (2021) COVID-19 PBMC atlas. It has
#' been preprocessed to include log-normalized counts, an author-provided
#' cell-type column, and PCA results, restricted to 5 shared cell types
#' (including CD14 monocytes) and a gene panel that always includes the
#' Yoshida et al. interferon-response signature.
#'
#' @details This dataset underwent the following steps:
#' \itemize{
#'  \item Downloads the processed reference (healthy) object from Zenodo
#'  (\url{https://doi.org/10.5281/zenodo.18274942}), itself derived from the
#'  Stephenson et al. (2021) COVID-19 PBMC atlas (CZI CELLxGENE).
#'  \item Restricts to 5 shared cell types (CD14 monocytes, CD4 T, CD8 T, B
#'  cells, NK_16hi) using the \code{author_cell_type_merged} column.
#'  \item Downsamples to at most 180 cells per cell type.
#'  \item Selects highly variable genes using the function \code{getTopHVGs}
#'  from the \code{scran} package, intersected between reference and query,
#'  always including the 25-gene Yoshida et al. (2019) interferon-response
#'  signature.
#'  \item Performs Principal Component Analysis (PCA) using the function
#'  \code{runPCA} from the \code{scater} package.
#' }
#'
#' @seealso Use \code{data("covid_reference_data")} to load and access the
#' resulting reference dataset, and \code{data("covid_query_data")} for the
#' corresponding query (severe COVID-19) dataset.
#'
#' @source Stephenson E, et al. (2021), processed and hosted at
#' \url{https://doi.org/10.5281/zenodo.18274942}.
#' @references Stephenson E, et al. (2021). Single-cell multi-omics analysis
#' of the immune response in COVID-19. Nature Medicine 27:904-916.
#'
#' @examples
#' # Load and explore the COVID-19 reference dataset
#' data("covid_reference_data")
#'
#' @keywords internal
#'
"covid_reference_data"

#' @title COVID-19 PBMC Query (Severe COVID-19) Single-Cell RNA-Seq Dataset
#'
#' @description This dataset contains a downsampled query (severe COVID-19
#' donor) subset of the Stephenson et al. (2021) COVID-19 PBMC atlas, for use
#' alongside \code{covid_reference_data}. Cell types are annotated using
#' Azimuth reference mapping (\code{azimuth_celltype_l1_merged}).
#'
#' @details See \code{covid_reference_data} for the shared processing steps.
#' CD14 monocytes, the focus of the accompanying case study vignette, are
#' downsampled to at most 450 cells (other cell types to at most 220) so
#' the interferon-activated subset described in the manuscript stays
#' well-represented.
#'
#' @seealso Use \code{data("covid_query_data")} to load and access the
#' resulting query dataset, and \code{data("covid_reference_data")} for
#' comparison with the reference dataset.
#'
#' @source Stephenson E, et al. (2021), processed and hosted at
#' \url{https://doi.org/10.5281/zenodo.18274942}.
#' @references Stephenson E, et al. (2021). Single-cell multi-omics analysis
#' of the immune response in COVID-19. Nature Medicine 27:904-916.
#'
#' @examples
#' # Load and explore the COVID-19 query dataset
#' data("covid_query_data")
#'
#' @keywords internal
#'
"covid_query_data"

#' @title MERFISH Colitis Reference (Healthy) Spatial Dataset
#'
#' @description This dataset contains a downsampled reference (Day 0,
#' healthy colon) subset of the Cadinu et al. (2024) mouse colitis MERFISH
#' dataset. It is a \linkS4class{SpatialExperiment} on a 943-gene targeted
#' panel, preprocessed to include log-normalized counts, a cell-type column,
#' and PCA results, restricted to 5 shared cell types (including
#' fibroblasts).
#'
#' @details This dataset underwent the following steps:
#' \itemize{
#'  \item Downloads the processed reference (Day 0) object from Zenodo
#'  (\url{https://doi.org/10.5281/zenodo.18274942}), itself derived from the
#'  Cadinu et al. (2024) MERFISH mouse colitis dataset (\code{MerfishData}
#'  package).
#'  \item Collapses inflammation-associated variants of a cell state (e.g.
#'  "Inflamed Fibroblast", present only in the query) into their parent
#'  lineage for a shared \code{cell_type_merged} column, while retaining
#'  the original fine-grained \code{tier2_merged} label.
#'  \item Restricts to 5 shared cell types (Fibroblast, Smooth Muscle,
#'  Epithelial, Other Immune, Endothelial).
#'  \item Downsamples to at most 220 cells per cell type.
#'  \item Performs Principal Component Analysis (PCA) using the function
#'  \code{runPCA} from the \code{scater} package on the full 943-gene panel.
#' }
#'
#' @seealso Use \code{data("merfish_reference_data")} to load and access the
#' resulting reference dataset, and \code{data("merfish_query_data")} for
#' the corresponding query (Day 9 colitis) dataset.
#'
#' @source Cadinu CA, et al. (2024), processed and hosted at
#' \url{https://doi.org/10.5281/zenodo.18274942}.
#' @references Cadinu CA, et al. (2024). Charting the cellular biogeography
#' in colitis reveals fibroblast trajectories and coordinated spatial
#' remodeling. Cell 187(8):2010-2028.
#'
#' @examples
#' # Load and explore the MERFISH reference dataset
#' data("merfish_reference_data")
#'
#' @keywords internal
#'
"merfish_reference_data"

#' @title MERFISH Colitis Query (Day 9 Colitis) Spatial Dataset
#'
#' @description This dataset contains a downsampled query (Day 9,
#' DSS-induced colitis at peak inflammation) subset of the Cadinu et al.
#' (2024) mouse colitis MERFISH dataset, for use alongside
#' \code{merfish_reference_data}.
#'
#' @details See \code{merfish_reference_data} for the shared processing
#' steps. Cells are downsampled to at most 260 per cell type. The
#' \code{tier2_merged} column retains the original fine-grained label (e.g.
#' distinguishing "Fibroblast" from "Inflamed Fibroblast"), which is not
#' passed to scDiagnostics functions directly but is useful for checking
#' which cells a diagnostic actually flags.
#'
#' @seealso Use \code{data("merfish_query_data")} to load and access the
#' resulting query dataset, and \code{data("merfish_reference_data")} for
#' comparison with the reference dataset.
#'
#' @source Cadinu CA, et al. (2024), processed and hosted at
#' \url{https://doi.org/10.5281/zenodo.18274942}.
#' @references Cadinu CA, et al. (2024). Charting the cellular biogeography
#' in colitis reveals fibroblast trajectories and coordinated spatial
#' remodeling. Cell 187(8):2010-2028.
#'
#' @examples
#' # Load and explore the MERFISH query dataset
#' data("merfish_query_data")
#'
#' @keywords internal
#'
"merfish_query_data"
