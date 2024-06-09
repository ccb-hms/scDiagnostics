#' @title Plot Projected Data on Discriminant Spaces
#'
#' @description 
#' This function plots the projected reference and query data on discriminant spaces.
#'
#' @details 
#' This function generates either a scatterplot or a boxplot to visualize the projected data onto the discriminant spaces.
#' For scatterplot, each point represents a projected data point, and colors are used to differentiate between different cell types 
#' and datasets. For boxplot, the distribution of the projected data values for each cell type is shown, separated by datasets.
#'
#' @param x An object of class \code{calculateDiscriminantSpace} containing the projected data on the discriminant space.. 
#' Each element of the list represents a combination of cell types and datasets. Each element should contain 'ref_proj' and 'query_proj' data frames.
#' @param cell_types A character vector specifying the cell types to plot. If not provided, all cell types will be plotted.
#' @param plot_type Type of plot to generate. Options are "scatterplot" and "boxplot". Default is "scatterplot".
#' @param ... Additional arguments to be passed to the plotting functions.
#'
#' @return A ggplot object representing the scatterplot or boxplot of the projected data.
#' 
#' @export
#' 
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#' 
#' @seealso \code{\link{calculateDiscriminantSpace}}
#' 
#' @rdname calculateDiscriminantSpace
#' 
# Function to plot the projected reference/query data on the discriminant spaces.
plot.calculateDiscriminantSpace <- function(x, cell_types, plot_type = c("scatterplot", "boxplot"), ...){
    
    # Check if query data is available in the object
    if(!("query_proj" %in% names(x[[names(x)[[1]]]])))
        stop("There is no query data to plot.")
    
    # Check input for plot_type
    plot_type <- match.arg(plot_type)
    if(!(plot_type %in% c("scatterplot", "boxplot"))){
        stop("The \'plot_type\' specified is not available.")
    }
    
    # Scatter plot
    if(plot_type == "scatterplot"){
        
        # Formatting of plot data
        plot_data <- rbind(data.frame(x[[1]][["ref_proj"]], data_type = rep("Reference", nrow(x[[1]][["ref_proj"]]))),
                           data.frame(x[[1]][["query_proj"]], data_type = rep("Query", nrow(x[[1]][["query_proj"]]))))
        # Create a new variable representing the combination of cell type and dataset
        plot_data[["cell_type_dataset"]] <- paste(plot_data[["data_type"]], plot_data[["cell_type"]], sep = " ")
        plot_data[["cell_type_combination"]] <- paste0(sort(unlist(strsplit(names(x)[1], "-"))), collapse = "-")
        full_data <- plot_data
        for(comb_id in 2:length(names(x))){
            # Formatting of plot data
            plot_data <- rbind(data.frame(x[[comb_id]][["ref_proj"]], data_type = rep("Reference", nrow(x[[comb_id]][["ref_proj"]]))),
                               data.frame(x[[comb_id]][["query_proj"]], data_type = rep("Query", nrow(x[[comb_id]][["query_proj"]]))))
            # Create a new variable representing the combination of cell type and dataset
            plot_data[["cell_type_dataset"]] <- paste(plot_data[["data_type"]], plot_data[["cell_type"]], sep = " ")
            plot_data[["cell_type_combination"]] <- paste0(sort(unlist(strsplit(names(x)[comb_id], "-"))), collapse = "-")
            full_data <- rbind(full_data, plot_data)
        }
        
        # Cell types 
        cell_types <- unique(unlist(unique(strsplit(names(x), "-"))))
        
        # Define the order of cell type and dataset combinations
        order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)),
                                    rep(sort(cell_types), each = 2))
        
        # Define the colors for cell types
        color_mapping <- setNames(RColorBrewer::brewer.pal(length(order_combinations), "Paired"), order_combinations)
        cell_type_colors <- color_mapping[order_combinations]
        
        # Reorder the levels of cell type and dataset factor
        full_data[["cell_type_dataset"]] <- factor(full_data[["cell_type_dataset"]], levels = order_combinations)

        # Generate scatter plot
        scatter_plot <- ggplot2::ggplot(full_data, ggplot2::aes(x = DV1, y = DV2, color = cell_type_dataset)) +
            ggplot2::geom_point(alpha = 0.5, size = 1) +
            ggplot2::scale_color_manual(values = cell_type_colors, name = "Cell Types") + 
            ggplot2::facet_wrap(~ cell_type_combination, scales = "free") +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                           plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                           axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                           axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10))
        return(scatter_plot)
        
    } else if (plot_type == "boxplot"){ # Boxplot
        
        # Check input for cell_types
        if(!(cell_types[1] %in% names(x))){
            
            cell_types <- paste0(rev(unlist(strsplit(cell_types[1], "-"))), collapse = "-")
            if(!(cell_types[1] %in% names(x)))
                stop("The specified combination \'cell_types\' is not available.")
        }
        
        # Formatting of plot data
        plot_data <- rbind(data.frame(x[[cell_types]]$ref_proj, data_type = rep("Reference", nrow(x[[cell_types]]$ref_proj))),
                           data.frame(x[[cell_types]]$query_proj, data_type = rep("Query", nrow(x[[cell_types]]$query_proj))))
        
        # Manually reshape the data to long format
        data_long <- data.frame(cell_type = rep(plot_data[["cell_type"]], times = 2),
                                data_type = rep(plot_data[["data_type"]], times = 2),
                                variable = rep(c("DV1", "DV2"), each = nrow(plot_data)),
                                value = c(plot_data$DV1, plot_data$DV2))
        
        # Create a new variable representing the combination of cell type and dataset
        data_long[["cell_type_dataset"]] <- paste(data_long[["data_type"]], data_long[["cell_type"]], sep = " ")
        
        # Cell types 
        cell_types <- unique(unlist(unique(strsplit(names(x), "-"))))
        
        # Define the order of cell type and dataset combinations
        order_combinations <- paste(rep(c("Reference", "Query"), length(cell_types)),
                                    rep(sort(cell_types), each = 2))
        
        # Reorder the levels of cell type and dataset factor
        data_long[["cell_type_dataset"]] <- factor(data_long[["cell_type_dataset"]], levels = order_combinations)

        # Define the colors for cell types
        color_mapping <- setNames(RColorBrewer::brewer.pal(length(order_combinations), "Paired"), order_combinations)
        cell_type_colors <- color_mapping[order_combinations]
        
        # Generate the boxplot
        box_plot <- ggplot2::ggplot(data_long, aes(x = cell_type, y = value, fill = cell_type_dataset)) +
            ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) + 
            ggplot2::facet_wrap(~ variable, scales = "free") +
            ggplot2::scale_fill_manual(values = cell_type_colors, name = "Cell Types") + 
            ggplot2::labs(x = "", y = "Value") +  
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_line(color = "gray", linetype = "dotted"),
                           plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                           axis.title = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 10),
                           axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10))
        return(box_plot)
    }
}