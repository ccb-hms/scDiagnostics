# -------------------------------------------
# Creation of Datasets - Multiple Cell Types
# -------------------------------------------

# Load data
sce <- scRNAseq::HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)

# Remove cells with NA cell type
sce <- sce[, !is.na(sce$reclustered.broad)]

# Divide the data into reference and query datasets
set.seed(100)
indices <- sample(ncol(SummarizedExperiment::assay(sce)),
                  size = floor(0.8 * ncol(SummarizedExperiment::assay(sce))),
                  replace = FALSE)
reference_data <- sce[, sample(indices, 1500)]
query_data <- sce[, -indices]

# log transform datasets
reference_data <- scuttle::logNormCounts(reference_data)
query_data <- scuttle::logNormCounts(query_data)

# Select specific column (cell) data
SummarizedExperiment::colData(reference_data) <-
    SummarizedExperiment::colData(reference_data)[, c("reclustered.broad"),
                                                  drop = FALSE]
SummarizedExperiment::colData(query_data) <-
    SummarizedExperiment::colData(query_data)[, c("percent.mito",
                                                  "reclustered.broad")]
names(SummarizedExperiment::colData(reference_data))[1] <- "expert_annotation"
names(SummarizedExperiment::colData(query_data))[1] <- "percent_mito"
names(SummarizedExperiment::colData(query_data))[2]  <- "expert_annotation"

# Selecting highly variable genes (can be customized by the user)
ref_var <- scran::getTopHVGs(reference_data, n = 500)
query_var <- scran::getTopHVGs(query_data, n = 500)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_var, query_var)
reference_data <- reference_data[common_genes, ]
query_data <- query_data[common_genes, ]

# Compute AUC gene set scores using AUCell with B cell signature
expression_matrix <- SummarizedExperiment::assay(query_data, "logcounts")
cells_rankings <- AUCell::AUCell_buildRankings(expression_matrix, plotStats = FALSE)
cd8_t_cell_signature <- c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK",
                          "GZMM", "PRF1", "NKG7", "CCL5", "CST7", "CTSW")

gene_sets <- list(CD8_T_cell_signature = cd8_t_cell_signature)
cells_AUC <- AUCell::AUCell_calcAUC(gene_sets, cells_rankings)
SummarizedExperiment::colData(query_data)$gene_set_scores <-
    SummarizedExperiment::assay(cells_AUC)["CD8_T_cell_signature", ]

# Run dimension reduction on the reference data
reference_data <- scater::runPCA(reference_data, ncomponents = 25)
reference_data <- scater::runTSNE(reference_data)
reference_data <- scater::runUMAP(reference_data)

# Run dimension reduction on the query data
query_data <- scater::runPCA(query_data, ncomponents = 25)
query_data <- scater::runTSNE(query_data)
query_data <- scater::runUMAP(query_data)

# Get cell type scores using SingleR (or any other cell type annotation method)
scores <- SingleR::SingleR(query_data, reference_data,
                           labels = reference_data$expert_annotation)
SummarizedExperiment::colData(query_data)$SingleR_annotation <- scores$labels
SummarizedExperiment::colData(query_data)$annotation_scores <- apply(scores$scores, 1, max)

# Remove counts assays
SummarizedExperiment::assays(reference_data) <-
    SummarizedExperiment::assays(reference_data)[-which(
        names(SummarizedExperiment::assays(reference_data)) == "counts")]
SummarizedExperiment::assays(query_data) <-
    SummarizedExperiment::assays(query_data)[-which(
        names(SummarizedExperiment::assays(query_data)) == "counts")]

# Save datasets to data/ folder
usethis::use_data(reference_data, compress = "xz", overwrite = TRUE)
usethis::use_data(query_data, compress = "xz", overwrite = TRUE)

