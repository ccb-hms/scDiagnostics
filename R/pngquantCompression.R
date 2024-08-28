#' @title Compress and Save a ggplot Object Using pngquant
#'
#' @description
#' This function saves a ggplot object to a file and compresses the resulting image using `pngquant`.
#'
#' @details
#' The function utilizes `ggsave` to save a `ggplot2` plot to a specified file path. After saving, it runs the `pngquant` command-line tool to compress the image. The compression quality can be controlled via the `quality` parameter. Ensure `pngquant` is installed and properly configured on your system for the compression to work.
#'
#' @param plot A `ggplot` object to be saved and compressed.
#' @param filename A character string specifying the file path where the plot will be saved.
#' @param height Numeric value indicating the height of the saved plot in inches. Default is 5 inches.
#' @param width Numeric value indicating the width of the saved plot in inches. Default is 10 inches.
#' @param quality A character string specifying the quality range for compression. The range should be in the form `"min-max"` where `min` and `max` are integers between 0 and 100. Default is `"1-5"`.
#' @param verbose Logical; if `TRUE`, additional information about the process is printed. Default is `FALSE`.
#'
#' @keywords internal
#'
#' @return The function does not return a value. It saves and compresses the plot directly to the specified file path.
#'
#' @author Anthony Christidis, \email{anthony-alexander_christidis@hms.harvard.edu}
#'
# Function to save and compress figure with pngquant
pngquantCompression <- function(plot, filename, height = 5, width = 10, quality = "1-5", verbose = FALSE) {
  
  # Save the plot using ggsave
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height)
  
  # Construct the system command for pngquant
  compress_cmd <- sprintf("pngquant --quality=%s --ext .png --force %s", 
                          quality, filename)
  
  # Execute the pngquant command to compress the image
  result <- system(compress_cmd, intern = TRUE, ignore.stderr = TRUE)
  
  # Check if compression with pngquant was executed
  if (length(result) == 0) {
    if(verbose){
      message("Plot saved and compressed at: ", filename)
    }
  } else {
    warning("Compression might have failed. Please check if pngquant is installed and properly configured.")
  }
}