# -----------------------------------------------
# Creation of Datasets - Zeisel Brain Benchmarking
# -----------------------------------------------

# Load data
sce <- scRNAseq::ZeiselBrainData()

# Log transform dataset
sce <- scuttle::logNormCounts(sce)

# Rename cell type column for clarity
names(SummarizedExperiment::colData(sce))[
    names(SummarizedExperiment::colData(sce)) == "level1class"] <- "true_cell_type"

# Divide the data into reference and query datasets (70/30 split)
set.seed(100)
indices <- sample(ncol(sce), size = floor(0.7 * ncol(sce)), replace = FALSE)
zeisel_reference_data <- sce[, indices]
zeisel_query_data <- sce[, -indices]

# Select specific column (cell) data
SummarizedExperiment::colData(zeisel_reference_data) <-
    SummarizedExperiment::colData(zeisel_reference_data)[, c("true_cell_type"),
                                                          drop = FALSE]
SummarizedExperiment::colData(zeisel_query_data) <-
    SummarizedExperiment::colData(zeisel_query_data)[, c("true_cell_type"),
                                                      drop = FALSE]

# Selecting highly variable genes (can be customized by the user)
ref_var <- scran::getTopHVGs(zeisel_reference_data, n = 250)
query_var <- scran::getTopHVGs(zeisel_query_data, n = 250)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_var, query_var)
zeisel_reference_data <- zeisel_reference_data[common_genes, ]
zeisel_query_data <- zeisel_query_data[common_genes, ]

# Run PCA on the reference and query data
zeisel_reference_data <- scater::runPCA(zeisel_reference_data, ncomponents = 10)
zeisel_query_data <- scater::runPCA(zeisel_query_data, ncomponents = 10)

# Remove counts assays
SummarizedExperiment::assays(zeisel_reference_data) <-
    SummarizedExperiment::assays(zeisel_reference_data)[-which(
        names(SummarizedExperiment::assays(zeisel_reference_data)) == "counts")]
SummarizedExperiment::assays(zeisel_query_data) <-
    SummarizedExperiment::assays(zeisel_query_data)[-which(
        names(SummarizedExperiment::assays(zeisel_query_data)) == "counts")]

# Normalize internal S4 representation to the current Bioconductor classes
zeisel_reference_data <- SingleCellExperiment::updateObject(zeisel_reference_data, check = FALSE)
zeisel_query_data <- SingleCellExperiment::updateObject(zeisel_query_data, check = FALSE)

# Save datasets to data/ folder
usethis::use_data(zeisel_reference_data, compress = "xz", overwrite = TRUE)
usethis::use_data(zeisel_query_data, compress = "xz", overwrite = TRUE)
