#' @title Plot Regression Results on Principal Components
#'
#' @description 
#' The S3 plot method generates plots to visualize the results of regression analyses 
#' performed on principal components (PCs) against cell types or dataset origin (query vs. reference).
#'
#' @details 
#' The S3 plot method generates, depending on the specified plot type, either the R-squared 
#' values or p-values resulting from the regression of principal components onto 
#' cell types or dataset origin (query vs. reference). For cell type regression, the plots show how well each 
#' PC correlates with different cell types. For dataset regression, the plots 
#' compare the PCs between query and reference datasets.
#'
#' @param x An object of class \code{regressPC} containing the output of the \code{regressPC} function
#' @param plot_type Type of plot to generate. Options are "r_squared" and "p-value". Default is "r-squared".
#' @param alpha Significance threshold p-values of coefficients. Default is 0.05.
#' @param ... Additional arguments to be passed to the plotting functions.
#'
#' @return The S3 plot method returns a \code{ggplot} object representing the specified plot type.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{regressPC}}
#' 
#' @rdname regressPC
#' 
# Function to plot the projected reference/query data on the discriminant spaces.
plot.regressPC <- function(x, plot_type = c("r_squared", "p-value"), alpha = 0.05, ...){
    
    # Check plot_type input
    plot_type <- match.arg(plot_type)
    if(!(plot_type %in% c("r_squared", "p-value")))
        stop("\'plot_type\' should be one of \'r_squared\' or \'p-value.\'")
    
    # Plots when regressing against cell types
    if(x[["indep_var"]] == "cell_type"){
        
        # R-squared plot
        if(plot_type == "r_squared"){
            plot_data <- data.frame(PC = paste0("PC", seq_len(length(x[["regression_summaries"]]))),
                                    r_squared = x[["r_squared"]])
            plot_data[["PC"]] <- factor(plot_data[["PC"]], levels = plot_data[["PC"]])
            
            plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]], y = .data[["r_squared"]], group = 1)) +
                ggplot2::geom_line() + ggplot2::geom_point() +
                ggplot2::labs(title = bquote(R^2 ~ of ~ "PC ~ " ~ "Cell Types"), x = NULL, y = bquote(R^2)) +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"))
        }
        
        # p-value plot
        if(plot_type == "p-value"){
            
            plot_data <- data.frame(PC = rep(names(x[["regression_summaries"]]), 
                                             each = nrow(x[["regression_summaries"]][[1]][["coefficients"]])), 
                                    cell_type = rownames(x[["regression_summaries"]][[1]][["coefficients"]]), p_value = NA)
            plot_data[["PC"]] <- factor(plot_data[["PC"]], levels = unique(plot_data[["PC"]]))
            
            # Define the order of cell type and dataset combinations
            cell_type_colors <- generateColors(sort(unique(plot_data[["cell_type"]])), paired = FALSE)
            
            for(pc in unique(plot_data[["PC"]])){
                plot_data[plot_data[["PC"]] == pc, "p_value"] <- x[["regression_summaries"]][[pc]][["coefficients"]][, "p.value"]
            }
            
            plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]], y = -log10(.data[["p_value"]]), 
                                                                color = .data[["cell_type"]], group = .data[["cell_type"]])) +
                ggplot2::geom_line() + ggplot2::geom_point() +
                ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Types") + 
                ggplot2::scale_y_continuous(trans = "log10") +
                ggplot2::labs(title = "P-values by Principal Component and Cell Type",
                              x = NULL, y = "-log10(p-value) of Cell Type") +
                ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"))
        }
        
    } else if (x[["indep_var"]] == "dataset"){
        
        # R-squared plot
        if(plot_type == "r_squared"){
            
            plot_data <- data.frame(PC = rep(names(x[[1]]), each = length(x) - 1), 
                                    cell_type = names(x)[-length(x)], r_squared = NA)
            plot_data[["PC"]] <- factor(plot_data[["PC"]], levels = unique(plot_data[["PC"]]))
            
            # Define the order of cell type and dataset combinations
            cell_type_colors <- generateColors(sort(unique(plot_data[["cell_type"]])), paired = FALSE)
            
            for(cell_type in unique(plot_data[["cell_type"]])){
                plot_data[plot_data[["cell_type"]] == cell_type, "r_squared"] <- 
                    unlist(lapply(seq_len(length(unique(plot_data[["PC"]]))), function(i) x[[cell_type]][[paste0("PC", i)]][["r_squared"]]))
            }
            
            plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]], y = .data[["r_squared"]], 
                                                                color = .data[["cell_type"]], group = .data[["cell_type"]])) +
                ggplot2::geom_line() + ggplot2::geom_point() +
                ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Type") + 
                ggplot2::labs(title = bquote(R^2 ~ of ~ "PC ~ " ~ "Dataset (Query Vs. Reference)"), x = NULL, y = bquote(R^2)) +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"))
        }
        
        # p-value plot
        if(plot_type == "p-value"){
            
            plot_data <- data.frame(PC = rep(names(x[[1]]), each = length(x) - 1), 
                                    cell_type = names(x)[-length(x)], p_value = NA)
            plot_data[["PC"]] <- factor(plot_data[["PC"]], levels = unique(plot_data[["PC"]]))
            
            # Define the order of cell type and dataset combinations
            cell_type_colors <- generateColors(sort(unique(plot_data[["cell_type"]])), paired = FALSE)
            
            for(pc in unique(plot_data[["PC"]])){
                plot_data[plot_data[["PC"]] == pc, "p_value"] <- unlist(lapply(unique(plot_data[["cell_type"]]), 
                                                                               function(t) x[[t]][[pc]][["coefficients"]][2, "p.value"]))
            }
            
            plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC"]], y = -log10(.data[["p_value"]]), 
                                                                color = .data[["cell_type"]], group = .data[["cell_type"]])) +
                ggplot2::geom_line() + ggplot2::geom_point() +
                ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Types") + 
                ggplot2::scale_y_continuous(trans = "log10") +
                ggplot2::labs(title = "P-values (Query Data Indicator) by Principal Component and Cell Type",
                              x = NULL, y = "-log10(p-value) of Query Indicator") +
                ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"))
        }
    }
    
    # Return plot object
    return(plot_obj)
}