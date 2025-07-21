#' @title Plot Principal Components for Different Cell Types
#'
#' @description This function generates a \code{ggplot2} visualization of principal components (PCs) for different
#' cell types across two datasets (query and reference), using either boxplots or violin plots.
#'
#' @details
#' The function \code{boxplotPCA} is designed to provide a visualization of principal component analysis (PCA) results. It projects
#' the query dataset onto the principal components obtained from the reference dataset. The results are then visualized
#' as boxplots or violin plots, grouped by cell types and datasets (query and reference). This allows for a comparative analysis of the
#' distributions of the principal components across different cell types and datasets. The function internally calls \code{projectPCA}
#' to perform the PCA projection. It then reshapes the output data into a long format suitable for ggplot2 plotting.
#'
#' @param query_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for the
#'                   query cells.
#' @param reference_data A \code{\linkS4class{SingleCellExperiment}} object containing numeric expression matrix for
#'                       the reference cells.
#' @param query_cell_type_col The column name in the \code{colData} of \code{query_data} that identifies the cell types.
#' @param ref_cell_type_col The column name in the \code{colData} of \code{reference_data} that identifies the cell types.
#' @param cell_types A character vector specifying the cell types to include in the plot. If NULL, all cell types are
#'                   included.
#' @param pc_subset A numeric vector specifying which principal components to include in the plot. Default is PC1 to PC5.
#' @param shape Character string indicating the plot type: "box" for boxplots or "violin" for violin plots.
#'              Default is "box".
#' @param assay_name Name of the assay on which to perform computations. Default is "logcounts".
#' @param max_cells Maximum number of cells to retain. If the object has fewer cells, it is returned unchanged.
#'                  Default is 2500.
#'
#' @return A ggplot object representing the boxplots or violin plots of specified principal components for the given
#'         cell types and datasets.
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
#' # Plot the PC data with boxplots (default)
#' pc_plot <- boxplotPCA(query_data = query_data,
#'                       reference_data = reference_data,
#'                       cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                       query_cell_type_col = "SingleR_annotation",
#'                       ref_cell_type_col = "expert_annotation",
#'                       pc_subset = 1:6)
#' pc_plot
#'
#' # Plot the PC data with violin plots
#' pc_violin <- boxplotPCA(query_data = query_data,
#'                         reference_data = reference_data,
#'                         cell_types = c("CD4", "CD8", "B_and_plasma", "Myeloid"),
#'                         query_cell_type_col = "SingleR_annotation",
#'                         ref_cell_type_col = "expert_annotation",
#'                         pc_subset = 1:6,
#'                         shape = "violin")
#' pc_violin
#'
#' @importFrom stats approxfun cancor density setNames
#' @importFrom utils combn
#'
# Function to plot PC for different cell types
boxplotPCA <- function(query_data,
                       reference_data,
                       query_cell_type_col,
                       ref_cell_type_col,
                       cell_types = NULL,
                       pc_subset = 1:5,
                       shape = c("box", "violin"),
                       assay_name = "logcounts",
                       max_cells = 2500){

    # Match the shape argument
    shape <- match.arg(shape)

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
                             query_cell_type_col = query_cell_type_col,
                             ref_cell_type_col = ref_cell_type_col,
                             pc_subset = pc_subset,
                             assay_name = assay_name)

    # Create the long format data frame manually
    pca_output <- pca_output[!is.na(pca_output[["cell_type"]]),]
    if(!is.null(cell_types)){
        if(all(cell_types %in% pca_output[["cell_type"]])){
            pca_output <- pca_output[which(pca_output[["cell_type"]] %in%
                                               cell_types),]
        } else{
            stop("One or more of the specified \'cell_types\' are not available.")
        }
    }
    pca_long <- data.frame(PC = rep(paste0("pc", pc_subset),
                                    each = nrow(pca_output)),
                           Value = unlist(c(pca_output[, pc_subset])),
                           dataset = rep(pca_output[["dataset"]],
                                         length(pc_subset)),
                           cell_type = rep(pca_output[["cell_type"]],
                                           length(pc_subset)))
    pca_long[["PC"]] <- toupper(pca_long[["PC"]])

    # Create a new variable representing the combination of cell type and dataset
    pca_long[["cell_type_dataset"]] <- paste(pca_long[["dataset"]],
                                             pca_long[["cell_type"]],
                                             sep = " ")

    # Define the order of cell type and dataset combinations
    order_combinations <- paste(rep(c("Reference", "Query"),
                                    length(unique(pca_long[["cell_type"]]))),
                                rep(sort(unique(pca_long[["cell_type"]])),
                                    each = 2))

    # Reorder the levels of cell type and dataset factor
    pca_long[["cell_type_dataset"]] <- factor(pca_long[["cell_type_dataset"]],
                                              levels = order_combinations)

    # Define the colors for cell types
    cell_type_colors <- generateColors(order_combinations,
                                       paired = TRUE)

    # Create the ggplot
    plot <- ggplot2::ggplot(pca_long, ggplot2::aes(
        x = .data[["cell_type"]],
        y = .data[["Value"]],
        fill = .data[["cell_type_dataset"]]))

    # Add either boxplot or violin plot based on the shape parameter
    if(shape == "box") {
        plot <- plot + ggplot2::geom_boxplot(alpha = 0.7,
                                             outlier.shape = NA,
                                             width = 0.7)
    } else { # shape == "violin"
        plot <- plot + ggplot2::geom_violin(alpha = 0.7,
                                            trim = FALSE,
                                            width = 0.7)
    }

    # Continue with common plot elements
    plot <- plot +
        ggplot2::facet_wrap(~ .data[["PC"]], scales = "free") +
        ggplot2::scale_fill_manual(values = cell_type_colors,
                                   name = "Cell Types") +
        ggplot2::labs(x = "", y = "PCA Score") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "white", color = "black", linewidth = 0.5),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "gray",
                                                     linetype = "dotted"),
            plot.title = ggplot2::element_text(size = 14,
                                               face = "bold", hjust = 0.5),
            axis.title = ggplot2::element_text(size = 12),
            axis.text = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(angle = 45,
                                                hjust = 1,
                                                size = 10))

    # Return the plot
    return(plot)
}
