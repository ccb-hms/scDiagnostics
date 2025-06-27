# ------------------------------
# Creation of Dataset - QC Data
# ------------------------------

# load reference dataset
reference_data <- celldex::fetchReference("hpca", "2024-02-26")

# Load query dataset (Bunis haematopoietic stem and progenitor cell data) from
# Bunis DG et al. (2021). Single-Cell Mapping of Progressive Fetal-to-Adult
# Transition in Human Naive T Cells Cell Rep. 34(1): 108573
qc_data <- scRNAseq::BunisHSPCData()
rownames(qc_data) <- rowData(qc_data)$Symbol

# Sample to reduce size of dataset
set.seed(100)
qc_data <- qc_data[, sample(seq_len(ncol(qc_data)), 750)]

# Only keep genes in reference data
common_genes <- intersect(rownames(qc_data),
                          rownames(reference_data))
qc_data <- qc_data[intersect(rownames(qc_data),
                             rownames(reference_data)),]
reference_data <- reference_data[intersect(rownames(qc_data),
                                           rownames(reference_data)),]

# Add QC metrics to query data
qc_data <- scuttle::addPerCellQCMetrics(qc_data)

# Log transform QC dataset
qc_data <- scuttle::logNormCounts(qc_data)

# Selecting highly variable genes (can be customized by the user)
qc_var <- scran::getTopHVGs(qc_data, n = 500)
qc_data <- qc_data[qc_var, ]

# Run SingleR to predict cell types
scores <- SingleR::SingleR(qc_data, reference_data, labels = reference_data$label.main)

# Assign predicted labels to query data
SummarizedExperiment::colData(qc_data)$SingleR_annotation <- scores$labels

# Assign scores to query data
SummarizedExperiment::colData(qc_data)$annotation_scores <- apply(scores$scores, 1, max)


# Select some cell data for final data
SummarizedExperiment::colData(qc_data) <-
    SummarizedExperiment::colData(qc_data)[, c("total",
                                               "SingleR_annotation",
                                               "annotation_scores")]

# Remove counts assays
SummarizedExperiment::assays(qc_data) <-
    SummarizedExperiment::assays(qc_data)[-which(
        names(SummarizedExperiment::assays(qc_data)) == "counts")]

# Save dataset to data/ folder
usethis::use_data(qc_data, compress = "xz", overwrite = TRUE)

