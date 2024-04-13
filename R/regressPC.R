#' Principal component regression
#'
#' This function performs linear regression of a covariate of interest
#' onto one or more principal components, based on the data in a
#' SingleCellExperiment object.
#'
#' @details Principal component regression, derived from PCA, can be
#'     used to quantify the variance explained by a covariate
#'     interest. Applications for single-cell analysis include
#'     quantification of batch removal, assessing clustering
#'     homogeneity, and evaluation of alignment of query and reference
#'     datasets in cell type annotation settings.  Briefly, the R^2 is
#'     calculated from a linear regression of the covariate B of
#'     interest onto each principal component. The variance
#'     contribution of the covariate effect per principal component is
#'     then calculated as the product of the variance explained by the
#'     ith principal component (PC) and the corresponding R2(PCi|B).
#'     The sum across all variance contributions by the covariate
#'     effects in all principal components gives the total variance
#'     explained by the covariate as follows:
#'
#'     Var(C|B) = sum_{i=1}^G Var(C|PC_i) * R^2 (PC_i | B)
#'
#'     where, Var(C|PCi) is the variance of the data matrix C explained
#'     by the ith principal component. See references.
#'
#' @param sce An object of class
#'     \code{\linkS4class{SingleCellExperiment}} containing the data
#'     for regression analysis.
#'
#' @param dep.vars character. Dependent variable(s). Determines which
#'     principal component(s) (e.g., "PC1", "PC2", etc.) are used as
#'     explanatory variables.  Principal components are expected to be
#'     stored in a PC matrix named \code{"PCA"} in the
#'     \code{reducedDims} of \code{sce}. Defaults to \code{NULL} which
#'     will then regress on each principal component present in the PC
#'     matrix.
#'
#' @param indep.var character. Independent variable. A column name in
#'     the \code{colData} of \code{sce} specifying the response
#'     variable.
#'
#' @param plot logical indicating whether to print two diagnostic
#'     plots: one showing the embedding of the two PCs with the
#'     highest R2, colored by the independent variable and a second
#'     showing the PC ~ independent variable R2 values.
#'
#' @param ... other arguments to pass to \link[scater]{plotPCA} for
#'     the second diagnostic plot.
#'
#' @return A \code{list} containing \itemize{ \item summaries of the
#'     linear regression models for each specified principal
#'     component, \item the corresponding R-squared (R2) values, \item
#'     the variance contributions for each principal component, and
#'     \item the total variance explained.}
#'
#' @references Luecken et al. Benchmarking atlas-level data
#'     integration in single-cell genomics. Nature Methods, 19:41-50,
#'     2022.
#'
#' @examples
#' library(scater)
#' library(scran)
#' library(scRNAseq)
#' library(SingleR)
#'
#' # Load data
#' sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)
#'
#' # Divide the data into reference and query datasets
#' set.seed(100)
#' indices <- sample(ncol(sce),
#'     size = floor(0.7 * ncol(sce)),
#'     replace = FALSE
#' )
#' ref <- sce[, indices]
#' query <- sce[, -indices]
#'
#' # log transform datasets
#' ref <- logNormCounts(ref)
#' query <- logNormCounts(query)
#'
#' # Run PCA
#' query <- runPCA(query)
#'
#' # Get cell type scores using SingleR
#' # Note: replace when using cell type annotation scores from other methods
#' scores <- SingleR(query, ref, labels = ref$reclustered.broad)
#'
#' # Add labels to query object
#' query$labels <- scores$labels
#'
#' # Specify the dependent variables (principal components) and
#' # independent variable (e.g., "labels")
#' dep.vars <- paste0("PC", 1:3)
#' indep.var <- "labels"
#'
#' # Perform linear regression on multiple principal components
#' res <- regressPC(
#'     sce = query,
#'     dep.vars = dep.vars,
#'     indep.var = indep.var
#' )
#'
#' # Obtain linear regression summaries and R-squared values
#' res$regression.summaries
#' res$rsquared
#'
#' @importFrom stats lm
#' @importFrom ggplot2 ggplot geom_point geom_line theme_bw ylim labs
#'     scale_x_continuous theme element_blank
#' @importFrom scater plotPCA
#' @import SingleCellExperiment
#' @export
#' 
regressPC <- function(sce,
                      dep.vars = NULL,
                      indep.var,
                      plot = TRUE,
                      ...) {
    ## sanity checks
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot("PCA" %in% reducedDimNames(sce))

    if (!is.null(dep.vars)) {
        stopifnot(all(dep.vars %in% colnames(reducedDim(sce, "PCA"))))
    }

    stopifnot(indep.var %in% colnames(colData(sce)))
    
    ## regress against all PCs if not instructed otherwise
    if (is.null(dep.vars))
        dep.vars <- colnames(reducedDim(sce, "PCA"))
    
    ## create a data frame with the dependent and independent variables
    df <- data.frame(
        Independent = sce[[indep.var]],
        reducedDim(sce, "PCA")[, dep.vars]
    )

    ## perform linear regression for each principal component
    .regress <- function(pc) {
        f <- paste0(pc, " ~ Independent")
        model <- lm(f, data = df)
        s <- summary(model)
        return(s)
    }
    summaries <- lapply(dep.vars, .regress)
    names(summaries) <- dep.vars

    ## calculate R-squared values
    rsq <- vapply(summaries, `[[`, numeric(1), x = "r.squared")

    ## calculate variance contributions by principal component
    ind <- match(dep.vars, colnames(reducedDim(sce, "PCA")))
    var.expl <- attr(reducedDim(sce, "PCA"), "percentVar")[ind]
    var.contr <- var.expl * rsq
    
    ## calculate total variance explained by summing the variance contributions
    total.var.expl <- sum(var.contr)

    ## return the summaries of the linear regression models,
    ## R-squared values, and variance contributions
    res <- list(
        regression.summaries = summaries,
        rsquared = rsq,
        var.contributions = var.contr,
        total.variance.explained = total.var.expl
    )

    if (plot) {
        message("Printing two diagnostic plots.")

        p1 <- plotPCA(
            sce,
            color_by = indep.var,
            ncomponents = base::rev(tail(order(rsq), 2)),
            ...
        )
        
        p2_input <- data.frame(
            x = dep.vars,
            i = seq_along(dep.vars),
            r2 = rsq
        )

        p2 <- ggplot(p2_input, aes(i, r2)) +
            geom_point() +
            geom_line() +
            theme_bw() +
            ylim(c(0, 1)) +
            labs(
                y = bquote(R^2 ~ of ~ "PC ~ " ~ .(indep.var))
            ) +
            scale_x_continuous(
                breaks = p2_input$i,
                labels = p2_input$x
            ) +
            theme(
                axis.title.x = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()
            )
        
        print(p1)
        print(p2)
    }

    res
}
