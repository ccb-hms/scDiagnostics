#' @title Perform Hotelling's T-squared Test on PCA Scores for Single-cell RNA-seq Data
#'
#' @description
#' Computes Hotelling's T-squared test statistic and p-values for each specified cell type
#' based on PCA-projected data from query and reference datasets.
#'
#' @details
#' This function calculates Hotelling's T-squared statistic for comparing multivariate means
#' between reference and query datasets, projected onto a subset of principal components (PCs).
#' It performs a permutation test to obtain p-values for each cell type specified.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the reference cells.
#' @param query_cell_type_col character. The column name in the \code{colData} of \code{query_data}
#' that identifies the cell types.
#' @param ref_cell_type_col character. The column name in the \code{colData} of \code{reference_data}
#' that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#' @param n_permutation Number of permutations to perform for p-value calculation. Default is 500.
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A named numeric vector of p-values from Hotelling's T-squared test for each cell type.
#'
#' @references
#' Hotelling, H. (1931). "The generalization of Student's ratio". *Annals of Mathematical Statistics*. 2 (3): 360–378. doi:10.1214/aoms/1177732979.
#'
#' @export
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
#' @examples
#' # Load data
#' data("reference_data")
#' data("query_data")
#'
#' # Get the p-values
#' p_values <- calculateHotellingPValue(query_data = query_data,
#'                                      reference_data = reference_data,
#'                                      query_cell_type_col = "SingleR_annotation",
#'                                      ref_cell_type_col = "expert_annotation",
#'                                      pc_subset = 1:10)
#' round(p_values, 5)
#'
# Function to perform Hotelling T^2 test for each cell type
# The test is performed on the PCA space of the reference data The query data projected onto PCA space of reference
calculateHotellingPValue <- function(query_data,
                                     reference_data,
                                     query_cell_type_col,
                                     ref_cell_type_col,
                                     cell_types = NULL,
                                     pc_subset = 1:5,
                                     n_permutation = 500,
                                     assay_name = "logcounts",
                                     max_cells = 2500) {

    # Check standard input arguments
    argumentCheck(query_data = query_data,
                  reference_data = reference_data,
                  query_cell_type_col = query_cell_type_col,
                  ref_cell_type_col = ref_cell_type_col,
                  cell_types = cell_types,
                  pc_subset_ref = pc_subset,
                  assay_name = assay_name)

    # Downsample query and reference data
    query_data <- downsampleSCE(sce = query_data,
                                max_cells = max_cells)
    reference_data <- downsampleSCE(sce = reference_data,
                                    max_cells = max_cells)

    # Get common cell types if they are not specified by user
    if(is.null(cell_types)){
        cell_types <- na.omit(unique(c(reference_data[[ref_cell_type_col]],
                                       query_data[[query_cell_type_col]])))
    }

    # Get the projected PCA data
    pca_output <- projectPCA(query_data = query_data,
                             reference_data = reference_data,
                             pc_subset = pc_subset,
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             assay_name = assay_name)
    cell_list <- split(pca_output, pca_output[["cell_type"]])

    # Perform permutation test with p-values
    p_values <- numeric(length(cell_types))
    names(p_values) <- cell_types
    for (cell_type in cell_types) {

        ref_ind <- cell_list[[cell_type]][["dataset"]] == "Reference"
        observed_t2_data <- as.matrix(cell_list[[cell_type]][, paste0("PC", pc_subset)])
        observed_t2 <- hotellingT2(observed_t2_data[ref_ind,],
                                   observed_t2_data[!ref_ind,])

        perm_t2 <- numeric(n_permutation)
        for(perm_id in seq_len(n_permutation)){

            ref_sample_id <- sample(seq_len(nrow(cell_list[[cell_type]])), sum(ref_ind), replace = FALSE)
            perm_t2_data <- as.matrix(cell_list[[cell_type]][, paste0("PC", pc_subset)])
            perm_t2[perm_id] <- hotellingT2(perm_t2_data[ref_sample_id,],
                                            perm_t2_data[-ref_sample_id,])
        }
        p_values[cell_type] <- mean(observed_t2 < perm_t2)
    }

    # Return p-values
    return(p_values)
}

#' @title Calculate Hotelling's T^2 Statistic
#'
#' @description
#' Calculates the Hotelling's T^2 statistic for comparing means of multivariate data.
#'
#' @param sample1 A numeric matrix or data frame of multivariate observations for sample 1, where rows are observations and columns
#' are variables.
#' @param sample2 A numeric matrix or data frame of multivariate observations for sample 2, with the same structure as sample 1.
#'
#' @keywords internal
#'
#' @return The Hotelling's T^2 statistic.
#'
# Function to calculate hotelling T2 statistic between two groups
hotellingT2 <- function(sample1, sample2) {

    # Number of observations in each sample
    n1 <- nrow(sample1)
    n2 <- nrow(sample2)

    # Number of variables (columns) in each sample
    p <- ncol(sample1)

    # Compute mean vectors for each sample
    mean1 <- colMeans(sample1)
    mean2 <- colMeans(sample2)

    # Compute covariance matrices for each sample
    cov1 <- cov(sample1)
    cov2 <- cov(sample2)

    # Pooled covariance matrix
    pooled_cov <- ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)

    # Compute Hotelling's T^2 statistic
    # T^2 = n1 * n2 / (n1 + n2) * (mean1 - mean2)' * pooled_cov^-1 * (mean1 - mean2)
    t2 <- n1 * n2 / (n1 + n2) * t(mean1 - mean2) %*% solve(pooled_cov) %*%
        (mean1 - mean2)

    # Return the computed T^2 statistic
    return(as.numeric(t2))
}

