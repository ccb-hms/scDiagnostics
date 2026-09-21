# ------------------------------------------------
# Creation of Datasets - COVID-19 PBMC Case Study
# ------------------------------------------------
# Source: Stephenson et al. (2021) COVID-19 PBMC atlas, processed and hosted
# at https://doi.org/10.5281/zenodo.18274942 (companion data for Christidis
# et al., Briefings in Bioinformatics 2026). Reference = healthy donors,
# query = severe COVID-19 donors, both already log-normalized with
# per-cell-type annotations from multiple tools (author-provided ground
# truth plus Azimuth/SingleR/CellTypist/scVI predictions).

library(SingleCellExperiment)

zenodo_base <- "https://zenodo.org/records/18274942/files/"
bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

normal_data <- readRDS(BiocFileCache::bfcrpath(bfc, paste0(zenodo_base, "normal_data_sce.rds")))
covid_data <- readRDS(BiocFileCache::bfcrpath(bfc, paste0(zenodo_base, "covid_data_sce.rds")))

# These objects were serialized by the manuscript's own (older) pipeline, so
# normalize their internal S4 representation to the current Bioconductor
# classes now, rather than relying on a lazy, version-dependent
# updateObject() the first time a (possibly newer) Bioconductor loads them
normal_data <- updateObject(normal_data, check = FALSE)
covid_data <- updateObject(covid_data, check = FALSE)

# Restrict to a handful of shared cell types for PCA context, matched by name
# between the author ground-truth (reference) and Azimuth (query) merged
# annotation columns
shared_cell_types <- c("CD14 mono", "CD4 T", "CD8 T", "B cell", "NK_16hi")

normal_data <- normal_data[, normal_data$author_cell_type_merged %in% shared_cell_types]
covid_data <- covid_data[, covid_data$azimuth_celltype_l1_merged %in% shared_cell_types]

# Downsample cells per type (CD14 monocytes, the focus of the analysis, get a
# higher cap in the query so the IFN-activated subset stays visible)
set.seed(100)
downsample_by_type <- function(sce, type_col, caps) {
    idx <- unlist(lapply(names(caps), function(ct) {
        ct_idx <- which(SummarizedExperiment::colData(sce)[[type_col]] == ct)
        sample(ct_idx, size = min(length(ct_idx), caps[[ct]]))
    }))
    sce[, idx]
}

ref_caps <- setNames(rep(180, length(shared_cell_types)), shared_cell_types)
query_caps <- setNames(rep(220, length(shared_cell_types)), shared_cell_types)
query_caps["CD14 mono"] <- 450

normal_data <- downsample_by_type(normal_data, "author_cell_type_merged", ref_caps)
covid_data <- downsample_by_type(covid_data, "azimuth_celltype_l1_merged", query_caps)

# Select specific column (cell) data
SummarizedExperiment::colData(normal_data) <-
    SummarizedExperiment::colData(normal_data)[, c("author_cell_type_merged"), drop = FALSE]
SummarizedExperiment::colData(covid_data) <-
    SummarizedExperiment::colData(covid_data)[, c("azimuth_celltype_l1_merged"), drop = FALSE]

# Selecting highly variable genes, forcing in the Yoshida et al. (2019)
# interferon-response signature used in the manuscript's gene-shift analysis
yoshida_ifn_signature <- c(
    "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L",
    "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2",
    "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22",
    "UBE2L6", "XAF1", "IRF7")

ref_var <- scran::getTopHVGs(normal_data, n = 500)
query_var <- scran::getTopHVGs(covid_data, n = 500)
common_genes <- union(intersect(ref_var, query_var), yoshida_ifn_signature)

covid_reference_data <- normal_data[common_genes, ]
covid_query_data <- covid_data[common_genes, ]

# Drop any pre-existing reduced dimensions (e.g. UMAP_scVI computed on the
# full dataset) so only a PCA computed on this downsampled panel remains
SingleCellExperiment::reducedDims(covid_reference_data) <- list()
SingleCellExperiment::reducedDims(covid_query_data) <- list()

# Run PCA on the reference and query data
covid_reference_data <- scater::runPCA(covid_reference_data, ncomponents = 15)
covid_query_data <- scater::runPCA(covid_query_data, ncomponents = 15)

# Remove counts assays
SummarizedExperiment::assays(covid_reference_data) <-
    SummarizedExperiment::assays(covid_reference_data)[-which(
        names(SummarizedExperiment::assays(covid_reference_data)) == "counts")]
SummarizedExperiment::assays(covid_query_data) <-
    SummarizedExperiment::assays(covid_query_data)[-which(
        names(SummarizedExperiment::assays(covid_query_data)) == "counts")]

# Save datasets to data/ folder
usethis::use_data(covid_reference_data, compress = "xz", overwrite = TRUE)
usethis::use_data(covid_query_data, compress = "xz", overwrite = TRUE)
