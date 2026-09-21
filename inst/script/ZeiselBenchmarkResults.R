# --------------------------------------------------------
# Creation of Dataset - Zeisel Anomaly Detection Benchmark
# --------------------------------------------------------
# Computes, once, the sensitivity/specificity/AUROC/AUPRC of detectAnomaly()
# and calculateReconstructionError() on the full Zeisel brain dataset across
# label-noise, class-imbalance, and batch-effect gradients, plus a
# hyperparameter grid search. Only the resulting (small) summary tables are
# bundled with the package; the full dataset used to compute them is not.

suppressMessages({
    library(scDiagnostics)
    library(SingleCellExperiment)
    library(scater)
    library(SingleR)
    library(scRNAseq)
    library(pROC)
    library(PRROC)
    library(dplyr)
})

# ______________________
# 1. DATA PREPARATION
# ______________________

sce <- scRNAseq::ZeiselBrainData()
sce <- logNormCounts(sce)
sce$true_cell_type <- sce$level1class

set.seed(0)
split_idx <- sample(seq_len(ncol(sce)), size = floor(ncol(sce) * 0.7))
base_ref <- sce[, split_idx]
base_query <- sce[, -split_idx]

# ___________________________
# 2. EVALUATION FUNCTIONS
# ___________________________

compute_metrics <- function(scores, truth_logical) {
    if (length(scores) == 0 || length(truth_logical) == 0) {
        return(c(AUROC = NA, AUPRC = NA))
    }
    valid_idx <- !is.na(scores) & !is.na(truth_logical)
    scores <- scores[valid_idx]
    truth_logical <- truth_logical[valid_idx]
    if (length(unique(truth_logical)) < 2) return(c(AUROC = NA, AUPRC = NA))
    roc_obj <- suppressMessages(roc(response = truth_logical, predictor = scores,
                                    quiet = TRUE))
    fg <- scores[truth_logical == TRUE]
    bg <- scores[truth_logical == FALSE]
    pr_obj <- suppressWarnings(pr.curve(scores.class0 = fg, scores.class1 = bg,
                                        curve = FALSE))
    c(AUROC = as.numeric(auc(roc_obj)), AUPRC = pr_obj$auc.integral)
}

calc_sens_spec <- function(preds, truth_logical) {
    TP <- sum(truth_logical == TRUE & preds == TRUE, na.rm = TRUE)
    TN <- sum(truth_logical == FALSE & preds == FALSE, na.rm = TRUE)
    FP <- sum(truth_logical == FALSE & preds == TRUE, na.rm = TRUE)
    FN <- sum(truth_logical == TRUE & preds == FALSE, na.rm = TRUE)
    c(Sensitivity = TP / (TP + FN), Specificity = TN / (TN + FP))
}

run_test <- function(ref_sce, query_sce, target_hidden_cluster) {
    pred <- SingleR(test = query_sce, ref = ref_sce, labels = ref_sce$true_cell_type)
    query_sce$SingleR_annotation <- pred$labels
    hidden_cells_predictions <-
        query_sce$SingleR_annotation[query_sce$true_cell_type == target_hidden_cluster]
    target_ref_cluster <- names(which.max(table(hidden_cells_predictions)))

    res_if <- suppressMessages(suppressWarnings(detectAnomaly(
        reference_data = ref_sce, query_data = query_sce,
        ref_cell_type_col = "true_cell_type", query_cell_type_col = "SingleR_annotation",
        cell_types = target_ref_cluster, pc_subset = NULL, n_hvgs = 100, n_tree = 500)))
    res_re <- suppressMessages(suppressWarnings(calculateReconstructionError(
        reference_data = ref_sce, query_data = query_sce,
        ref_cell_type_col = "true_cell_type", query_cell_type_col = "SingleR_annotation",
        cell_types = target_ref_cluster, pc_subset = 1:5, n_hvgs = 100)))

    query_target_cells <- query_sce[, query_sce$SingleR_annotation == target_ref_cluster]
    valid_cells <- query_target_cells$true_cell_type %in%
        c(target_ref_cluster, target_hidden_cluster)
    truth_logical <- query_target_cells$true_cell_type[valid_cells] == target_hidden_cluster

    if_scores <- res_if[[target_ref_cluster]]$query_anomaly_scores[valid_cells]
    re_scores <- res_re[[target_ref_cluster]]$query_reconstruction_errors[valid_cells]
    if_preds <- res_if[[target_ref_cluster]]$query_anomaly[valid_cells]
    re_preds <- res_re[[target_ref_cluster]]$query_anomaly[valid_cells]

    m_if <- compute_metrics(if_scores, truth_logical)
    m_re <- compute_metrics(re_scores, truth_logical)
    ss_if <- calc_sens_spec(if_preds, truth_logical)
    ss_re <- calc_sens_spec(re_preds, truth_logical)

    data.frame(
        Method = c("Isolation_Forest", "Reconstruction_Error"),
        AUROC = c(m_if["AUROC"], m_re["AUROC"]),
        AUPRC = c(m_if["AUPRC"], m_re["AUPRC"]),
        Sensitivity = c(ss_if["Sensitivity"], ss_re["Sensitivity"]),
        Specificity = c(ss_if["Specificity"], ss_re["Specificity"]))
}

# ______________________________________
# 3. BASELINE + GRADIENT BENCHMARK
# ______________________________________

all_results <- list()

message("Running baselines...")
res_dist <- run_test(base_ref[, base_ref$true_cell_type != "astrocytes_ependymal"],
                     base_query, "astrocytes_ependymal")
res_dist$Test <- "Distinct (Astrocytes)"; res_dist$TestGroup <- "Baseline"
all_results[[length(all_results) + 1]] <- res_dist

res_rel <- run_test(base_ref[, base_ref$true_cell_type != "pyramidal SS"],
                    base_query, "pyramidal SS")
res_rel$Test <- "Related (Pyramidal SS)"; res_rel$TestGroup <- "Baseline"
all_results[[length(all_results) + 1]] <- res_rel

res_rare <- run_test(base_ref[, base_ref$true_cell_type != "microglia"],
                     base_query, "microglia")
res_rare$Test <- "Rare (Microglia)"; res_rare$TestGroup <- "Baseline"
all_results[[length(all_results) + 1]] <- res_rare

run_gradient_suite <- function(hidden_cluster, group_name,
                                imb_levels = c(300, 100, 50, 20, 10)) {
    message(paste("Running gradients for:", group_name))
    ref_grad <- base_ref[, base_ref$true_cell_type != hidden_cluster]

    pred <- SingleR(test = base_query, ref = ref_grad, labels = ref_grad$true_cell_type)
    hidden_preds <- pred$labels[base_query$true_cell_type == hidden_cluster]
    target_ref <- names(which.max(table(hidden_preds)))

    # Noise
    for (nl in c(0.0, 0.1, 0.2, 0.3, 0.4)) {
        ref_tmp <- ref_grad
        n_shuffle <- floor(ncol(ref_tmp) * nl)
        if (n_shuffle > 0) {
            shuffle_idx <- sample(seq_len(ncol(ref_tmp)), size = n_shuffle)
            ref_tmp$true_cell_type[shuffle_idx] <- sample(ref_tmp$true_cell_type[shuffle_idx])
        }
        res_tmp <- run_test(ref_tmp, base_query, hidden_cluster)
        res_tmp$Test <- "Noise"; res_tmp$X_Value <- as.character(nl)
        res_tmp$TestGroup <- group_name
        all_results[[length(all_results) + 1]] <<- res_tmp
    }

    # Imbalance
    target_ref_idx <- which(ref_grad$true_cell_type == target_ref)
    max_clean_size <- floor(length(target_ref_idx) / 10) * 10
    valid_levels <- sort(unique(sapply(imb_levels, min, max_clean_size)), decreasing = TRUE)

    for (cl in valid_levels) {
        ref_tmp <- ref_grad
        ca1_idx <- which(ref_tmp$true_cell_type == target_ref)
        keep_ca1 <- sample(ca1_idx, size = cl)
        ref_tmp <- ref_tmp[, -setdiff(ca1_idx, keep_ca1)]

        res_tmp <- run_test(ref_tmp, base_query, hidden_cluster)
        res_tmp$Test <- "Imbalance"; res_tmp$X_Value <- as.character(cl)
        res_tmp$TestGroup <- group_name
        all_results[[length(all_results) + 1]] <<- res_tmp
    }

    # Batch effect
    for (bl in c(0.0, 0.5, 1.0, 1.5, 2.0)) {
        query_tmp <- base_query
        if (bl > 0) {
            genes_to_shift <- sample(seq_len(nrow(query_tmp)), size = floor(nrow(query_tmp) * 0.2))
            temp_counts <- as.matrix(logcounts(query_tmp))
            temp_counts[genes_to_shift, ] <- temp_counts[genes_to_shift, ] + bl
            logcounts(query_tmp) <- temp_counts
        }
        res_tmp <- run_test(ref_grad, query_tmp, hidden_cluster)
        res_tmp$Test <- "Batch"; res_tmp$X_Value <- as.character(bl)
        res_tmp$TestGroup <- group_name
        all_results[[length(all_results) + 1]] <<- res_tmp
    }
}

set.seed(1)
run_gradient_suite("pyramidal SS", "Related", imb_levels = c(300, 100, 50, 20, 10))
run_gradient_suite("microglia", "Rare", imb_levels = c(300, 100, 50, 20, 10))

zeisel_gradient_results <- bind_rows(all_results)

# ___________________________________
# 4. ISOLATION FOREST TUNING GRID
# ___________________________________

message("Running isolation forest tuning grid...")

ref_sce <- scDiagnostics::processPCA(base_ref, n_hvgs = 1000)
ref_sce <- ref_sce[, ref_sce$true_cell_type != "pyramidal SS"]

pred <- SingleR(test = base_query, ref = ref_sce, labels = ref_sce$true_cell_type)
query_sce <- base_query
query_sce$SingleR_annotation <- pred$labels

target_cluster <- "pyramidal CA1"
query_target_cells <- query_sce[, query_sce$SingleR_annotation == target_cluster]
valid_cells <- query_target_cells$true_cell_type %in% c("pyramidal CA1", "pyramidal SS")
is_true_anomaly <-
    query_target_cells$true_cell_type[valid_cells] == "pyramidal SS"

pc_list <- list("1:3" = 1:3, "1:5" = 1:5, "1:10" = 1:10)
hvg_list <- c(50, 100, 500)
thresholds <- list(
    list(name = "Abs (0.5)", method = "absolute", mad_val = 2, abs_val = 0.5),
    list(name = "MAD (2x)", method = "MAD", mad_val = 2, abs_val = 0.5),
    list(name = "MAD (3x)", method = "MAD", mad_val = 3, abs_val = 0.5))

calc_grid_metrics <- function(config_name, mode, param, thresh_name, preds) {
    ss <- calc_sens_spec(preds, is_true_anomaly)
    data.frame(Config = config_name, Mode = mode, Parameter = param,
              Threshold = thresh_name, Sensitivity = ss["Sensitivity"],
              Specificity = ss["Specificity"])
}

zeisel_if_tuning_results <- data.frame()
for (pc_name in names(pc_list)) {
    for (thresh in thresholds) {
        res <- suppressMessages(suppressWarnings(detectAnomaly(
            reference_data = ref_sce, query_data = query_sce,
            ref_cell_type_col = "true_cell_type", query_cell_type_col = "SingleR_annotation",
            cell_types = target_cluster, pc_subset = pc_list[[pc_name]], n_hvgs = 1000,
            n_tree = 300, threshold_method = thresh$method, mad_multiplier = thresh$mad_val,
            anomaly_threshold = thresh$abs_val)))
        preds <- res[[target_cluster]]$query_anomaly[valid_cells]
        zeisel_if_tuning_results <- rbind(zeisel_if_tuning_results,
            calc_grid_metrics(paste0("PCA_", pc_name), "Global PCA", pc_name, thresh$name, preds))
    }
}
for (hvg in hvg_list) {
    for (thresh in thresholds) {
        res <- suppressMessages(suppressWarnings(detectAnomaly(
            reference_data = ref_sce, query_data = query_sce,
            ref_cell_type_col = "true_cell_type", query_cell_type_col = "SingleR_annotation",
            cell_types = target_cluster, pc_subset = NULL, n_hvgs = hvg, n_tree = 300,
            threshold_method = thresh$method, mad_multiplier = thresh$mad_val,
            anomaly_threshold = thresh$abs_val)))
        preds <- res[[target_cluster]]$query_anomaly[valid_cells]
        zeisel_if_tuning_results <- rbind(zeisel_if_tuning_results,
            calc_grid_metrics(paste0("HVG_", hvg), "Targeted HVGs", as.character(hvg), thresh$name, preds))
    }
}
rownames(zeisel_if_tuning_results) <- NULL

# _____________________________________
# 5. RECONSTRUCTION ERROR TUNING GRID
# _____________________________________

message("Running reconstruction error tuning grid...")

pc_list_re <- list("1:3" = 1:3, "1:5" = 1:5, "1:10" = 1:10)
hvg_list_re <- c(50, 100, 500)
mad_list_re <- c(2, 3)

zeisel_re_tuning_results <- data.frame()
for (hvg in hvg_list_re) {
    for (pc_name in names(pc_list_re)) {
        pc <- pc_list_re[[pc_name]]
        if (max(pc) >= hvg) next
        for (mad_val in mad_list_re) {
            res <- suppressMessages(suppressWarnings(calculateReconstructionError(
                reference_data = ref_sce, query_data = query_sce,
                ref_cell_type_col = "true_cell_type", query_cell_type_col = "SingleR_annotation",
                cell_types = target_cluster, pc_subset = pc, n_hvgs = hvg,
                mad_multiplier = mad_val)))
            preds <- res[[target_cluster]]$query_anomaly[valid_cells]
            ss <- calc_sens_spec(preds, is_true_anomaly)
            zeisel_re_tuning_results <- rbind(zeisel_re_tuning_results, data.frame(
                HVGs = as.character(hvg), PCs = pc_name,
                MAD_Threshold = paste0("MAD (", mad_val, "x)"),
                Sensitivity = ss["Sensitivity"], Specificity = ss["Specificity"]))
        }
    }
}
zeisel_re_tuning_results$ConfigLabel <-
    paste0(zeisel_re_tuning_results$HVGs, " (", zeisel_re_tuning_results$PCs, ")")
zeisel_re_tuning_results$HVGs <-
    factor(zeisel_re_tuning_results$HVGs, levels = as.character(sort(hvg_list_re)))
zeisel_re_tuning_results$PCs <-
    factor(zeisel_re_tuning_results$PCs, levels = names(pc_list_re))
rownames(zeisel_re_tuning_results) <- NULL

# ______________
# 6. SAVE DATA
# ______________

zeisel_benchmark_results <- list(
    gradients = zeisel_gradient_results,
    if_tuning = zeisel_if_tuning_results,
    re_tuning = zeisel_re_tuning_results)

usethis::use_data(zeisel_benchmark_results, compress = "xz", overwrite = TRUE)
