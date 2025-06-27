#' @title Generate Paired Colors for Cell Types
#'
#' @description
#' This function assigns paired colors (light and dark) to a list of cell type names. The colors are selected from various color
#' palettes in the `pals` package.
#'
#' @details
#' The function uses color palettes from the `pals` package to generate colors or pairs of colors (light and dark) for each cell
#' type name provided. It cycles through different color families (blues, greens, reds, oranges, purples, purd and greys) to create
#' the colors
#'
#' @param cell_type_names A character vector of cell type names that need to be assigned colors.
#' @param paired If TRUE, the colored returned should be paired. Default is FALSE.
#'
#' @keywords internal
#'
#' @return A named character vector where the names are the original cell type names, and the values are the assigned colors.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to generate paired colors
generateColors <- function(cell_type_names, paired = FALSE){

    if (isTRUE(paired)) {

        if(length(cell_type_names) > 20)
            stop("There can be at most 10 cell types for plot function.")

        full_colors <- c(c("#BDD7E7", "#3182BD"),
                         c("#FCAE91", "#DE2D26"),
                         c("#BAE4B3", "#31A354"),
                         c("#FED766", "#E68E00"),
                         c("#CBC9E2", "#756BB1"),
                         c("#D7B5D8", "#DD1C77"),
                         c("#80e0d0", "#20b2aa"),
                         c("#D2B48C", "#8B4513"),
                         c("#FFD700", "#B8860B"),
                         c("#B7CE95", "#6B8E23"))[seq_len(
                             length(cell_type_names))]
        names(full_colors) <- cell_type_names

    } else if (isFALSE(paired)){

        if(length(cell_type_names) > 10)
            stop("There can be at most 10 cell types for plot function.")

        full_colors <- c("#74C476",
                         "#FB6A4A",
                         "#6BAED6",
                         "#FDB32C",
                         "#9E9AC8",
                         "#DF65B0",
                         "#50C9BD",
                         "#AF7D50",
                         "#FFC400",
                         "#A2B964")[seq_len(length(cell_type_names))]
        names(full_colors) <- cell_type_names
    }
    return(full_colors)
}



