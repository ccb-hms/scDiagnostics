scDiagnostics: diagnostic functions to assess the quality of cell type
assignments in single-cell RNA-seq
================
Smriti Chawla
2023-07-08

<h1>
Getting started
</h1>

# Installation

<h4>
Install a Bioconductor Package
</h4>

To install a package from Bioconductor use a following command:

``` r
BiocManager::install("scDiagnostics")
```

To build the package vignettes upon installation use:

``` r
BiocManager::install("ccb-hms/scDiagnostics",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

To install a package directly from a GitHub repository, you will need a
remotes package and command is as follows:

``` r
remotes::install_github("ccb-hms/scDiagnostics")
```

# Prerequisites

R version tested: 4.2.3 (2023-03-15)

To use the scDiagnostics package, you need the following R packages
installed:

- SingleCellExperiment
- scRNAseq
- scater
- scran
- ggplot2
- AUCell
- SingleR
- Matrix
- corrplot
- RColorBrewer

# Usage

To explore the capabilities of the scDiagnostics package, you can load
your own data or use the provided example with publicly available data
of Marrow tissue single cell gene expression profiles from He S et
al. (2020). Single-cell transcriptome profiling of an adult human cell
atlas of 15 major organs. Genome Biol 21, 1:294. This dataset is
obtained from the scRNAseq R package and is used for analyzing and
visualizing diagnostic plots to inspect query/test and reference
datasets and check the appropriateness of cell type assignments. In the
provided example, the dataset is divided into 70% reference and 30%
query data.

``` r
## Loading libraries
library(scDiagnostics)
library(scater)
library(scran)
library(scRNAseq)
library(RColorBrewer)
library(SingleR)
library(AUCell)
library(corrplot)
```

## Scatter Plot: QC stats vs. Annotation Scores

Let’s load the Marrow dataset for demonstration purposes. The dataset
represents single-cell gene expression profiles of Marrow tissue. In
this example, we will analyze the relationship between user-defined QC
stats and annotation scores for a specific cell type.

To perform this analysis, we first divide the dataset into reference and
query datasets. The reference dataset serves as a reference for cell
type annotation, while the query dataset contains cells that we want to
assign cell types to. We then log-transform the expression values in
both datasets for further analysis.

Next, we use the SingleR package to obtain cell type scores for the
query dataset. The scatter plot is generated to examine the relationship
between the percentage of mitochondrial genes and SingleR scores. It
allow users to visulaize relation of QC stats for all cell types or cell
types of interest.

It’s important to note that the use of SingleR for cell type annotation
is just one example. Users can use any other cell type annotation method
to obtain cell type scores. The scatter plot provides a visual
representation of this relationship, allowing for further examination
and interpretation of the cell type assignments in the dataset.

``` r

   # Load data
   sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)

   # Divide the data into reference and query datasets
   indices <- sample(ncol(assay(sce)), size = floor(0.7 * ncol(assay(sce))), replace = FALSE)
   ref_data <- sce[, indices]
   query_data <- sce[, -indices]

   # Log-transform datasets
   ref_data <- logNormCounts(ref_data)
   query_data <- logNormCounts(query_data)

   # Run PCA
   ref_data <- runPCA(ref_data)
   query_data <- runPCA(query_data)

   # Get cell type scores using SingleR
   pred <- SingleR(query_data, ref_data, labels = ref_data$reclustered.broad)
   pred <- as.data.frame(pred)
   
   # Assign labels to query data
   colData(query_data)$labels <- pred$labels
   
   # get scores
   scores <- apply(pred[,1:4], 1, max)

   # Assign scores to query data
   colData(query_data)$cell_scores <- scores

   # Generate scatter plots
   p1 <- plotQCvsAnnotation(query_data, "percent.mito", "labels", "cell_scores", c("CD4", "CD8"))
   p1 + xlab("percent.mito")
```

<img src="man/figures/Scatter plot QC stats vs annotation scores-1.png" width="100%" />

Scatter plot for visualizing relationship between percentage of
mitochondrial gene and cell annotation scores for the cell types.

``` r
   p2 <- plotQCvsAnnotation(query_data, "percent.mito", "labels", "cell_scores", NULL)
   p2 + xlab("percent.mito")
```

<img src="man/figures/Scatter plot QC stats vs annotation scores all cell types-1.png" width="100%" />

## Examining Distribution of QC stats and Annotation Scores

In addition to the scatter plot, we can gain further insights into the
gene expression profiles by visualizing the distribution of user defined
QC stats and annotation scores for all the cell types or specific cell
types. This allows us to examine the variation and patterns in
expression levels and scores across cells assigned to the cell type of
interest.

To accomplish this, we create two separate histograms. The first
histogram displays the distribution of the annotation scores. The x-axis
represents the scores, while the y-axis represents the frequency of
cells with a given score.

The second histogram visualizes the distribution of QC stats. This
provides insights into the overall gene expression levels for the
specific cell type. Here in this particular example we are investigating
percentage of mitochondrial genes.

By examining the histograms, we can observe the range, shape, and
potential outliers in the distribution of both annotation scores and
library size. This allows us to assess the appropriateness of the cell
type assignments and identify any potential discrepancies or patterns in
the gene expression profiles for the specific cell type.

``` r
# Generate histogram
p <- histQCvsAnnotation(query_data, "percent.mito", "labels", "cell_scores", NULL)

# Access the histogram of QC stats (the first element in the list)
qc_histogram <- p[[1]]
qc_histogram
```

<img src="man/figures/Distribution of library size and Annotation Score-1.png" width="100%" />

``` r

# Access the histogram of annotation scores (the second element in the list)
scores_histogram <- p[[2]]
scores_histogram
```

<img src="man/figures/Distribution of library size and Annotation Score-2.png" width="100%" />

The example code provided demonstrates how to utilize the
plotCellTypeDistribution function with the necessary data and packages.
Here we are visualizing the distribution of counts and scores for the
“CD4” cell type.

## Exploring Gene Expression Distribution

To gain insights into the gene expression values for a specific gene and
its distribution across cells, we can use the
plotGeneExpressionDistribution function. This function allows us to
visualize the distribution of expression values for a particular gene of
interest, both overall and within specific cell types.

``` r
# Generate histogram
plotGeneExpressionDistribution(query_data, "labels", "B_and_plasma", "VPREB3")
```

<img src="man/figures/histogram gene expression-1.png" width="100%" />

In the provided example, we are examining the distribution of expression
values for the gene “VPREB3” in the dataset. The function generates a
histogram to display the distribution of expression values across all
cells, with the x-axis representing the expression values and the y-axis
representing the frequency of cells within each expression range.
Additionally, the expression values can be further stratified by
specific cell types, providing insights into the gene expression
patterns within different cell populations.

## Visualize Gene Expression on Dimensional Reduction Plot

o gain insights into the gene expression patterns and their
representation in a dimensional reduction space, we can utilize the
plotGeneExpressionDimred function. This function allows us to plot the
gene expression values of a specific gene on a dimensional reduction
plot generated using methods like t-SNE, UMAP, or PCA. Each single cell
is color-coded based on its expression level of the gene of interest.

In the provided example, we are visualizing the gene expression values
of the gene “VPREB3” on a PCA plot. The PCA plot represents the cells in
a lower-dimensional space, where the x-axis corresponds to the first
principal component (Dimension 1) and the y-axis corresponds to the
second principal component (Dimension 2). Each cell is represented as a
point on the plot, and its color reflects the expression level of the
gene “VPREB3,” ranging from low (lighter color) to high (darker color).

``` r
# Generate dimension reduction plot color code by gene expression
plotGeneExpressionDimred(query_data, "PCA", c(1, 2), "VPREB3")
```

<img src="man/figures/scatter plot gene expression-1.png" width="100%" />

The dimensional reduction plot allows us to observe how the gene
expression of “VPREB3” is distributed across the cells and whether any
clusters or patterns emerge in the data. It provides a visual
representation of the gene expression landscape and can help identify
potential relationships or subpopulations based on the expression of the
gene of interest.

## Visualize Gene Sets or Pathway Scores on Dimensional Reduction Plot

In addition to examining individual gene expression patterns, it is
often useful to assess the collective activity of gene sets or pathways
within single cells. This can provide insights into the functional
states or biological processes associated with specific cell types or
conditions. To facilitate this analysis, the scDiagnostics package
includes a function called plotGeneSetScores that enables the
visualization of gene set or pathway scores on a dimensional reduction
plot.

The plotGeneSetScores function allows you to plot gene set or pathway
scores on a dimensional reduction plot generated using methods such as
PCA, t-SNE, or UMAP. Each single cell is color-coded based on its scores
for specific gene sets or pathways. This visualization helps identify
the heterogeneity and patterns of gene set or pathway activity within
the dataset, potentially revealing subpopulations with distinct
functional characteristics.

``` r

# Compute scores using AUCell
expression_matrix <- assay(query_data, "logcounts")
cells_rankings <- AUCell_buildRankings(expression_matrix, plotStats = F)

# Generate gene sets
gene_set1 <- sample(rownames(expression_matrix), 10)
gene_set2 <- sample(rownames(expression_matrix), 20)

gene_sets <- list(geneSet1 = gene_set1,
                  geneSet2 = gene_set2)

# Calculate AUC scores for gene sets
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

# Assign scores to colData
colData(query_data)$geneSetScores <- assay(cells_AUC)["geneSet1", ]

# Plot gene set scores on PCA
plotGeneSetScores(query_data, method = "PCA", feature = "geneSetScores")
```

<img src="man/figures/Visualize gene set or pathway scores on dimensional reduction scatter plot -1.png" width="100%" />

In the provided example, we demonstrate the usage of the
plotGeneSetScores function using the AUCell package to compute gene set
or pathway scores. Custom gene sets are generated for demonstration
purposes, but users can provide their own gene set scores using any
method of their choice. It is important to ensure that the scores are
assigned to the colData of the query_data object and specify the correct
feature name for visualization.

By visualizing gene set or pathway scores on a dimensional reduction
plot, you can gain a comprehensive understanding of the functional
landscape within your single-cell gene expression dataset and explore
the relationships between gene set activities and cellular phenotypes.

## Visualizing Reference and Query Cell Types using Multidimensional Scaling (MDS)

Here, we perform Multidimensional Scaling (MDS) analysis on the query
and reference datasets to examine their similarity. The dissimilarity
matrix is calculated based on the correlation between the datasets, and
MDS is used to obtain low-dimensional coordinates for each cell.
Subsequently, a scatter plot is generated, where each data point
represents a cell, and the cell types are color-coded using custom
colors provided by the user. This visualization enables the comparison
of cell type distributions between the query and reference datasets in a
reduced-dimensional space.

``` r
## Selcting highly variable genes
ref_var <- getTopHVGs(ref_data, n=2000)
query_var <- getTopHVGs(query_data, n=2000)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_var, query_var)

# Select desired cell types
selected_cell_types <- c("CD4", "CD8", "B_and_plasma")
ref_data_subset <- ref_data[common_genes, ref_data$reclustered.broad %in% selected_cell_types]
query_data_subset <- query_data[common_genes, query_data$reclustered.broad %in% selected_cell_types]

# Extract cell types for visualization
ref_labels <- ref_data_subset$reclustered.broad
query_labels <- query_data_subset$reclustered.broad

# Combine the cell type labels from both datasets
mdata <- c(paste("Query", query_labels), paste("Reference", ref_labels))

## Define the cell types and legend order
cell_types <- c("Query CD8", "Reference CD8", "Query CD4", "Reference CD4", "Query B_and_plasma", "Reference B_and_plasma")
legend_order <- cell_types

## Define the colors for cell types
color_palette <- brewer.pal(length(cell_types), "Paired")
color_mapping <- setNames(color_palette, cell_types)
cell_type_colors <- color_mapping[cell_types]

## Generate the MDS scatter plot with cell type coloring
visualizeCellTypeMDS(query_data_subset, ref_data_subset, mdata, cell_type_colors, legend_order)
```

<img src="man/figures/CMD scatter plot-1.png" width="100%" />

Upon examining the MDS scatter plot, we observe that the “CD4” and “CD8”
cell types overlap to some extent.By observing the proximity or overlap
of different cell types, we can gain insights into their potential
relationships or shared characteristics.

The selection of custom genes and desired cell types depends on the
user’s research interests and goals. It allows for flexibility in
focusing on specific genes and examining particular cell types of
interest in the visualization.

## Cell Type-specific Pairwise Correlation Analysis and Visualization

This analysis aims to explore the correlation patterns between different
cell types in a single-cell gene expression dataset. It involves
comparing the gene expression profiles of cells from a reference dataset
and a query dataset.

To perform the analysis, we start by computing the pairwise correlations
between the query and reference cells for selected cell types (“CD4”,
“CD8”, “B_and_plasma”). The Spearman correlation method is used, user
can also use Pearsons correlation coeefficient.

This will return average correlation matrix which can be visulaized by
user’s method of choice. Here, the results are visualized as a
correlation plot using the corrplot package.

``` r
selected_cell_types <- c("CD4", "CD8", "B_and_plasma")
cor_matrix_avg <- computeAveragePairwiseCorrelation(query_data_subset, ref_data_subset, "labels", "reclustered.broad", selected_cell_types, "spearman")

# Plot the pairwise average correlations using corrplot
corrplot(cor_matrix_avg, method = "number", tl.col = "black")
```

<img src="man/figures/Cell Type-specific Pairwise Correlation Analysis and Visualization -1.png" width="100%" />

This analysis allows us to examine the correlation patterns between
different cell types in the single-cell gene expression dataset. By
visualizing the average pairwise correlations, we can gain insights into
the relationships and similarities between cell types based on their
gene expression profiles.

In this case, users have the flexibility to extract the gene expression
profiles of specific cell types from the reference and query datasets
and provide these profiles as input to the function. Additionally, they
can select their own set of genes that they consider relevant for
computing the pairwise correlations. For demonstartion we have used
common highly variable genes from reference and query dataset.

By providing their own gene expression profiles and choosing specific
genes, users can focus the analysis on the cell types and genes of
interest to their research question.

## Pairwise Distance Analysis and Density Visualization

his function allows for the calculation of pairwise distances between
query and reference cells of a specific cell type in a single-cell gene
expression dataset. The pairwise distances are calculated using a
specified distance metric, such as “euclidean”, “manhattan”, etc.
Subsequently, the function generates density plots to visualize the
distribution of distances for different pairwise comparisons.

The analysis focuses on the specific cell type “CD8”, and the pairwise
distances are calculated using the “euclidean” distance metric.

``` r
calculatePairwiseDistancesAndPlotDensity(query_data_subset, ref_data_subset, "labels", "reclustered.broad", "CD8", "euclidean")
```

<img src="man/figures/Pairwise Distance Analysis and Density Visualization-1.png" width="100%" />

Further, user can also use correlation for calculation of pairwise
distances between query and reference cells of a specific cell type.

The analysis focuses on the specific cell type “CD8”, and the pairwise
distances are calculated using the correlation distance metric. User can
use spearman or pearson correlation coefficient as a method of choice.

``` r
calculatePairwiseDistancesAndPlotDensity(query_data_subset, ref_data_subset, "labels", "reclustered.broad", "CD8", "correlation" ,"spearman")
```

<img src="man/figures/Pairwise Distance Analysis and Density Visualization correlation based-1.png" width="100%" />

By utilizing this function, users can explore the pairwise distances
between query and reference cells of a specific cell type and gain
insights into the distribution of distances through density plots. This
analysis aids in understanding the similarities and differences in gene
expression profiles for the selected cell type within the query and
reference datasets.

## Linear regression analysis

Performing linear regression analysis on a SingleCellExperiment object
enables users to examine the relationship between a principal component
(PC) from the dimension reduction slot and an independent variable of
interest. By specifying the desired dependent variable as one of the
principal components (e.g., “PC1”, “PC2”, etc.) and providing the
corresponding independent variable from the colData slot, users can
explore the associations between these variables within the single-cell
gene expression dataset (reference and query).

``` r
summary <- performLinearRegression(query_data, "PC1", "labels")
#> 
#> Call:
#> lm(formula = Dependent ~ Independent, data = df)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -8.4684 -2.4142 -0.7017  2.0922 15.5537 
#> 
#> Coefficients:
#>                    Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)         -8.2479     0.2892  -28.52   <2e-16 ***
#> IndependentCD4       4.7333     0.3441   13.76   <2e-16 ***
#> IndependentCD8      13.8447     0.3351   41.32   <2e-16 ***
#> IndependentMyeloid   9.4121     0.5908   15.93   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 3.494 on 965 degrees of freedom
#> Multiple R-squared:  0.7014, Adjusted R-squared:  0.7005 
#> F-statistic: 755.7 on 3 and 965 DF,  p-value: < 2.2e-16
print(summary)
#> 
#> Call:
#> lm(formula = Dependent ~ Independent, data = df)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -8.4684 -2.4142 -0.7017  2.0922 15.5537 
#> 
#> Coefficients:
#>                    Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)         -8.2479     0.2892  -28.52   <2e-16 ***
#> IndependentCD4       4.7333     0.3441   13.76   <2e-16 ***
#> IndependentCD8      13.8447     0.3351   41.32   <2e-16 ***
#> IndependentMyeloid   9.4121     0.5908   15.93   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 3.494 on 965 degrees of freedom
#> Multiple R-squared:  0.7014, Adjusted R-squared:  0.7005 
#> F-statistic: 755.7 on 3 and 965 DF,  p-value: < 2.2e-16
```

By conducting linear regression, one can assess whether the PC values
are significantly associated with the cell types. This analysis helps
uncover whether there is a systematic variation in PC values across
different cell types. Understanding the relationship between PC values
and cell types can provide valuable insights into the biological or
technical factors driving cellular heterogeneity. It can help identify
PC dimensions that capture variation specific to certain cell types or
distinguish different cellular states.

# Conclusion

In this analysis, we have demonstrated the capabilities of the
scDiagnostics package for assessing the appropriateness of cell
assignments in single-cell gene expression profiles. By utilizing
various diagnostic functions and visualization techniques, we have
explored different aspects of the data, including total UMI counts,
annotation scores, gene expression distributions, dimensional reduction
plots, gene set scores, pairwise correlations, pairwise distances, and
linear regression analysis.

Through the scatter plots, histograms, and dimensional reduction plots,
we were able to gain insights into the relationships between gene
expression patterns, cell types, and the distribution of cells in a
reduced-dimensional space. The examination of gene expression
distributions, gene sets, and pathways allowed us to explore the
functional landscape and identify subpopulations with distinct
characteristics within the dataset. Additionally, the pairwise
correlation and distance analyses provided a deeper understanding of the
similarities and differences between cell types, highlighting potential
relationships and patterns.

------------------------------------------------------------------------

## R.session Info

    R version 4.3.1 (2023-06-16)
    Platform: aarch64-apple-darwin20 (64-bit)
    Running under: macOS Ventura 13.5

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] corrplot_0.92               AUCell_1.22.0              
     [3] SingleR_2.2.0               RColorBrewer_1.1-3         
     [5] scRNAseq_2.14.0             scran_1.28.1               
     [7] scater_1.29.2               ggplot2_3.4.2              
     [9] scuttle_1.9.4               scDiagnostics_0.99.0       
    [11] SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2
    [13] Biobase_2.60.0              GenomicRanges_1.52.0       
    [15] GenomeInfoDb_1.36.0         IRanges_2.34.0             
    [17] S4Vectors_0.38.1            BiocGenerics_0.46.0        
    [19] MatrixGenerics_1.12.0       matrixStats_1.0.0          

    loaded via a namespace (and not attached):
      [1] rstudioapi_0.15.0             magrittr_2.0.3               
      [3] ggbeeswarm_0.7.2              GenomicFeatures_1.52.1       
      [5] farver_2.1.1                  rmarkdown_2.22               
      [7] BiocIO_1.10.0                 zlibbioc_1.46.0              
      [9] vctrs_0.6.3                   Rsamtools_2.16.0             
     [11] memoise_2.0.1                 DelayedMatrixStats_1.22.0    
     [13] RCurl_1.98-1.12               htmltools_0.5.5              
     [15] S4Arrays_1.0.4                progress_1.2.2               
     [17] AnnotationHub_3.8.0           curl_5.0.1                   
     [19] BiocNeighbors_1.18.0          cachem_1.0.8                 
     [21] GenomicAlignments_1.36.0      igraph_1.5.0                 
     [23] mime_0.12                     lifecycle_1.0.3              
     [25] pkgconfig_2.0.3               rsvd_1.0.5                   
     [27] Matrix_1.6-0                  R6_2.5.1                     
     [29] fastmap_1.1.1                 GenomeInfoDbData_1.2.10      
     [31] shiny_1.7.4                   digest_0.6.31                
     [33] colorspace_2.1-0              AnnotationDbi_1.62.1         
     [35] dqrng_0.3.0                   irlba_2.3.5.1                
     [37] ExperimentHub_2.8.0           RSQLite_2.3.1                
     [39] beachmat_2.16.0               labeling_0.4.2               
     [41] filelock_1.0.2                fansi_1.0.4                  
     [43] httr_1.4.6                    compiler_4.3.1               
     [45] bit64_4.0.5                   withr_2.5.0                  
     [47] BiocParallel_1.34.2           viridis_0.6.3                
     [49] DBI_1.1.3                     highr_0.10                   
     [51] R.utils_2.12.2                biomaRt_2.56.1               
     [53] rappdirs_0.3.3                DelayedArray_0.26.3          
     [55] rjson_0.2.21                  bluster_1.10.0               
     [57] tools_4.3.1                   vipor_0.4.5                  
     [59] beeswarm_0.4.0                interactiveDisplayBase_1.38.0
     [61] httpuv_1.6.11                 R.oo_1.25.0                  
     [63] glue_1.6.2                    restfulr_0.0.15              
     [65] promises_1.2.0.1              grid_4.3.1                   
     [67] cluster_2.1.4                 generics_0.1.3               
     [69] gtable_0.3.3                  R.methodsS3_1.8.2            
     [71] ensembldb_2.24.0              data.table_1.14.8            
     [73] hms_1.1.3                     xml2_1.3.5                   
     [75] BiocSingular_1.16.0           ScaledMatrix_1.8.1           
     [77] metapod_1.8.0                 utf8_1.2.3                   
     [79] XVector_0.40.0                ggrepel_0.9.3                
     [81] BiocVersion_3.17.1            pillar_1.9.0                 
     [83] stringr_1.5.0                 limma_3.56.2                 
     [85] later_1.3.1                   dplyr_1.1.2                  
     [87] BiocFileCache_2.8.0           lattice_0.21-8               
     [89] rtracklayer_1.60.0            bit_4.0.5                    
     [91] annotate_1.78.0               tidyselect_1.2.0             
     [93] locfit_1.5-9.8                Biostrings_2.68.1            
     [95] knitr_1.43                    gridExtra_2.3                
     [97] ProtGenerics_1.32.0           edgeR_3.42.4                 
     [99] xfun_0.39                     statmod_1.5.0                
    [101] stringi_1.7.12                lazyeval_0.2.2               
    [103] yaml_2.3.7                    evaluate_0.21                
    [105] codetools_0.2-19              tibble_3.2.1                 
    [107] graph_1.78.0                  BiocManager_1.30.21          
    [109] cli_3.6.1                     xtable_1.8-4                 
    [111] munsell_0.5.0                 Rcpp_1.0.10                  
    [113] dbplyr_2.3.2                  png_0.1-8                    
    [115] XML_3.99-0.14                 parallel_4.3.1               
    [117] ellipsis_0.3.2                blob_1.2.4                   
    [119] prettyunits_1.1.1             AnnotationFilter_1.24.0      
    [121] sparseMatrixStats_1.12.0      bitops_1.0-7                 
    [123] GSEABase_1.62.0               viridisLite_0.4.2            
    [125] scales_1.2.1                  purrr_1.0.1                  
    [127] crayon_1.5.2                  rlang_1.1.1                  
    [129] cowplot_1.1.1                 KEGGREST_1.40.0              
