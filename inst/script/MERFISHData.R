# ---------------------------------------------------
# Creation of Datasets - MERFISH Spatial Case Study
# ---------------------------------------------------
# Source: Cadinu et al. (2024) mouse colitis MERFISH dataset, processed and
# hosted at https://doi.org/10.5281/zenodo.18274942 (companion data for
# Christidis et al., Briefings in Bioinformatics 2026). Reference = healthy
# colon (Day 0), query = DSS-induced colitis at peak inflammation (Day 9).
# Both are SpatialExperiment objects on the same 943-gene targeted panel.

library(SpatialExperiment)

zenodo_base <- "https://zenodo.org/records/18274942/files/"
bfc <- BiocFileCache::BiocFileCache(ask = FALSE)

healthy_data <- readRDS(BiocFileCache::bfcrpath(bfc, paste0(zenodo_base, "healthy_data.rds")))
dss9_data <- readRDS(BiocFileCache::bfcrpath(bfc, paste0(zenodo_base, "dss9_data.rds")))

# These objects were serialized by the manuscript's own (older) pipeline, so
# normalize their internal S4 representation to the current Bioconductor
# classes now, rather than relying on a lazy, version-dependent
# updateObject() the first time a (possibly newer) Bioconductor loads them
healthy_data <- updateObject(healthy_data, check = FALSE)
dss9_data <- updateObject(dss9_data, check = FALSE)

# The Day 9 (query) data further splits some cell states into an
# inflammation-associated variant (e.g. "Inflamed Fibroblast") that does not
# exist at Day 0 (reference). Collapse those back into their parent lineage
# for the shared cell-type column used by scDiagnostics, while keeping the
# original fine-grained label so the ground truth can be checked afterwards.
collapse_inflamed <- function(x) {
    x <- as.character(x)
    x[x == "Inflamed Fibroblast"] <- "Fibroblast"
    x[x == "Inflamed SMC"] <- "Smooth Muscle"
    x[x == "Inflamed Epithelial"] <- "Epithelial"
    x
}
healthy_data$cell_type_merged <- collapse_inflamed(healthy_data$tier2_merged)
dss9_data$cell_type_merged <- collapse_inflamed(dss9_data$tier2_merged)

# Restrict to a handful of shared cell types for context (Fibroblast is the
# focus of the analysis)
shared_cell_types <- c("Fibroblast", "Smooth Muscle", "Epithelial",
                       "Other Immune", "Endothelial")
healthy_data <- healthy_data[, healthy_data$cell_type_merged %in% shared_cell_types]
dss9_data <- dss9_data[, dss9_data$cell_type_merged %in% shared_cell_types]

# Downsample cells per type
set.seed(100)
downsample_by_type <- function(sce, type_col, max_n) {
    idx <- unlist(lapply(unique(SummarizedExperiment::colData(sce)[[type_col]]), function(ct) {
        ct_idx <- which(SummarizedExperiment::colData(sce)[[type_col]] == ct)
        sample(ct_idx, size = min(length(ct_idx), max_n))
    }))
    sce[, idx]
}

healthy_data <- downsample_by_type(healthy_data, "cell_type_merged", 220)
dss9_data <- downsample_by_type(dss9_data, "cell_type_merged", 260)

# Select specific column (cell) data (keep the fine-grained label too, for
# checking how many of the ground-truth "Inflamed Fibroblast" cells are
# actually flagged as anomalous)
SummarizedExperiment::colData(healthy_data) <-
    SummarizedExperiment::colData(healthy_data)[, c("cell_type_merged", "tier2_merged")]
SummarizedExperiment::colData(dss9_data) <-
    SummarizedExperiment::colData(dss9_data)[, c("cell_type_merged", "tier2_merged")]

# Drop any pre-existing reduced dimensions so only a freshly computed PCA
# remains
SingleCellExperiment::reducedDims(healthy_data) <- list()
SingleCellExperiment::reducedDims(dss9_data) <- list()

# Run PCA on the reference and query data (fixed 943-gene panel, no HVG
# reduction since it is already a targeted platform panel)
merfish_reference_data <- scater::runPCA(healthy_data, ncomponents = 15)
merfish_query_data <- scater::runPCA(dss9_data, ncomponents = 15)

# Remove counts assays
SummarizedExperiment::assays(merfish_reference_data) <-
    SummarizedExperiment::assays(merfish_reference_data)[-which(
        names(SummarizedExperiment::assays(merfish_reference_data)) == "counts")]
SummarizedExperiment::assays(merfish_query_data) <-
    SummarizedExperiment::assays(merfish_query_data)[-which(
        names(SummarizedExperiment::assays(merfish_query_data)) == "counts")]

# Save datasets to data/ folder
usethis::use_data(merfish_reference_data, compress = "xz", overwrite = TRUE)
usethis::use_data(merfish_query_data, compress = "xz", overwrite = TRUE)
