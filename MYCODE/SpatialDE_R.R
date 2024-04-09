# Mar. 18, 2024
# SpatialDE

# SpatialDE is a method to identify genes which significantly depend on spatial coordinates in
# non-linear and non-parametric ways. The intended applications are spatially resolved RNA-seq

# SpatialDE also provides automatic expression histology, a method that groups genes into common
# spatial patterns (and conversely reveal histological patterns based on gene expression).

# The key features of the method are:
# Unsupervised - No need to define spatial regions
# Non-parametric and non-linear expression patterns
# Automatic histology based on spatially co-expressed genes
# Extremely fast - Transcriptome wide tests takes only a few mins on normal computers

### install packages
# BiocManager::install("spatialDE")
# BiocManager::install("DropletUtils")

### libraries
library(spatialDE)
library(ggplot2)
library(SpatialExperiment)

### 1.1 Load data
# Example: Mouse Olfactory Bulb
data("Rep11_MOB_0")
data("MOB_sample_info")

# Rep11_MOB_0 contains spatial expression data for 16218 genes on 262 spots
dim(Rep11_MOB_0) # 16218   262
Rep11_MOB_0[1:5, 1:5]
# Column represents the coordinates (x-axis * y-axis)

# MOB_sample_info contains a dataframe with coordinates for each spot
dim(MOB_sample_info) # 260   3
head(MOB_sample_info)

# Filter out pratically unobserved genes
sum(rowSums(Rep11_MOB_0) >= 3) # 14859

Rep11_MOB_0 <- Rep11_MOB_0[rowSums(Rep11_MOB_0) >= 3, ]
dim(Rep11_MOB_0) # 14859 262
head(Rep11_MOB_0, 2) # gene by spot matrix

# Get total counts for every spots
diff <- setdiff(colnames(Rep11_MOB_0), rownames(MOB_sample_info))
diff

Rep11_MOB_0 <- Rep11_MOB_0[, row.names(MOB_sample_info)]
MOB_sample_info$total_counts <- colSums(Rep11_MOB_0)
head(MOB_sample_info)

# Get coordinates from MOB_sample_info
X <- MOB_sample_info[, c("x", "y")]
head(X)

### 1.2 Stabilize
# SpatialDE method assumes normally distributed data

# stabilize() is used to stabilize the variance of negative binomial distribution data
# and returns the stabilized count matrix.
normal_expression <- stabilize(Rep11_MOB_0)

range(normal_expression) # 0.8810934 6.2174312
normal_expression[1:5, 1:5]
Rep11_MOB_0[1:5, 1:5]

### 1.3 Regress out
# Regress out the effect of library size.

resid_expression <- regress_out(
  counts = normal_expression, # matrix of variance stabilized counts, resulting from stabilize()
  sample_info = MOB_sample_info # df with samples as rows and at least a column with total_counts
) # return a matrix of normalized counts

resid_expression[1:5, 1:5]

### 1.4 Run
# To reduce the running time, SpatialDE test is run on a subset of 1000 genes.
# Running it on the complete dataset takes about 10 mins.

# Perform SpatialDE test
sample_resid_expression <- head(resid_expression, 1000)

results <- spatialDE::run(
  sample_resid_expression, 
  coordinates = X,
  verbose = TRUE
)
# g: the name of the gene
# pval: the p-value for spatial differential expression
# qval: Significance after correcting for multiple testing
# l: A parameter indicating the distance scale a gene changes expression over

# run() is used to perform spatial differential expression testing and returns a dataframe
# containing the differential expression results of genes and related statistical information

head(results[order(results$qval), ])

### 1.5 Model search
# Classify the DE genes to interpetable DE classes
# Apply the model search on filted DE results
# using the threshold of 0.05 for the q-value
de_results <- results[results$qval < 0.05, ]
head(de_results)
dim(de_results)

ms_results <- model_search(
  # matrix-like object of normalized counts. E.g. resulting from regress_out()
  sample_resid_expression, 
  coordinates = X, 
  de_results = de_results # dataframe resulting from run()
)

head(ms_results)
dim(ms_results)

### 1.6 Spatial patterns
# Group spatial variable genes (SVGs) into spatial patterns
# using automatic expression histology (AEH)

# AEH: 一种自动化的组织学表达分析技术
# 结合了计算机视觉和机器学习技术，用于对组织学样本中的表达模式进行分析
# To perform AEH, the genes should be filtered by SpatialDE significance

sp <- spatial_patterns(
  sample_resid_expression,
  coordinates = X,
  de_results = de_results,
  qval_thresh = 0.05, # only rows in de_results with qval < qval_thresh will be kept
  n_patterns = 4L, # the number of spatial patterns (4)
  length = 1.5 # the characteristic length scale of the clusters
)

dim(sp$pattern_results)
# spatial_patterns() is used to group spatially variable genes into spatial patterns
# and returns pattern membership information and the posterior mean expression of genes in the pattern

# results: list of 2 dataframe (pattern_results, patterns)
# pattern_results: pattern membership information for each gene
# patterns: the posterior mean underlying expression for genes in given spatial patterns
sp$pattern_results

# order by pattern
sp$pattern_results[order(sp$pattern_results$pattern), ]

### 1.7 Plots

# Single gene (one of the most significant genes)
gene <- "Pcp4"

ggplot(data = MOB_sample_info, aes(x = x, y = y, color = normal_expression[gene, ])) +
  geom_point(size = 7) +
  ggtitle(gene) +
  scale_color_viridis_c() +
  labs(color = gene)

gene <- "X2900097C17Rik"

ggplot(data = MOB_sample_info, aes(x = x, y = y, color = normal_expression[gene, ])) +
  geom_point(size = 7) +
  ggtitle(gene) +
  scale_color_viridis_c() +
  labs(color = gene)

# Multiple genes
ordered_de_results <- de_results[order(de_results$qval), ]

multiGenePlots(normal_expression,
               coordinates = X,
               ordered_de_results[1:6, ]$g,
               point_size = 4,
               viridis_option = "D",
               dark_theme = FALSE
)

# Plot fraction spatial variance (FSV) vs q-value
# 空间方差比与Q值之间的关系图
# FSV是一种用于衡量在空间上的变异程度的指标
# FSV表示该基因表达变异中有多少是由于空间位置造成的
FSV_sig(results = results, ms_results = ms_results)
# A higher FSV value means that the expression pattern of the gene has significant variability in space
# indicating that the expression of the gene is greatly affected by the spatial location



### 2. Spatial Experiment Integration

# For SpatialExperiment object, we need to transpose the counts matrix
# in order to have genes on rows and spot on columns. 
# For this example, run spatialDE on the first 1000 genes

partial_counts <- head(Rep11_MOB_0, 1000)

# The SpatialExperiment class is designed to represent spatially resolved trancriptomics (ST) data
# supports storage of spatial information via "spatialCoords"
#          storage of images via "imgData"
# Extends the SingleCellExperiment class
spe <- SpatialExperiment(
  assays = list(counts = partial_counts), # assays: primary and transformed data
  spatialData = DataFrame(MOB_sample_info[, c("x", "y")]),
  spatialCoordsNames = c("x", "y") # default
)

out <- spatialDE(spe, assay_type = "counts", verbose = FALSE)
head(out[order(out$qval), ])

### 2.1 Plot spatial patterns of multiple genes (using SpatialExperiment object)
spe_results <- out[out$qval < 0.05, ]

ordered_spe_results <- spe_results[order(spe_results$qval), ]

multiGenePlots(
  spe,
  assay_type = "counts",
  ordered_spe_results[1:6, ]$g,
  point_size = 4,
  viridis_option = "D",
  dark_theme = FALSE
)

### 2.2 Classify spatially variable genes with model_search and spatial_patterns
msearch <- modelSearch(
  spe,
  de_results = out,
  qval_thresh = 0.05,
  verbose = FALSE
)
head(msearch)

spatterns <- spatialPatterns(
  spe,
  de_results = de_results,
  qval_thresh = 0.05,
  n_patterns = 4L,
  length = 1.5,
  verbose = FALSE
)

spatterns$pattern_results



#####################################
### SpatialExperiment

### libraries
library(DropletUtils)

# 1.1 Load data
example(read10xVisium, echo = FALSE)

# SpatialExperimental class
# contains imgData
spe

# spatialCoords are stored inside the "int_colData" 
# and are directly accessible via the corresponding accessor
head(spatialCoords(spe))

# The corresponding column names
spatialCoordsNames(spe)

# All image related data are stored inside the "int_metadata's" imgData field as a df
imgData(spe)

spi <- getImg(spe)
spi

identical(spi, imgData(spe)$data[[1]])
# [1] TRUE

# Data available in an object of class SpatialImage may be accessed 
# via the imgRaster() and imgSource() accessors
plot(imgRaster(spe))
imgSource(spe)

### Other operations
# Add or remove images using addImg() and rmvImg()
# Image transformations like rotation and mirroring







