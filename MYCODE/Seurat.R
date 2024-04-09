# Mar. 8, 2024

### install packages
install.packages("Seurat")
install.packages("languageserver")
# Install the remotes package
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes") # nolint
# }
# install.packages('Signac')
# remotes::install_github("satijalab/seurat-data", quiet = TRUE) # nolint
# remotes::install_github("satijalab/azimuth", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

### libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# PBMC: Peripheral Blood Mononuclear Cells 

### 1. Setup the Seurat object
pbmc.data <- Read10X(data.dir = "Data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
# specify minimum number of cells per gene and minimum number of features(genes) per cell
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# the size of the dense matrix
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

# the size of the sparse matrix
# Seurat typically uses sparse matrices to store single-cell sequencing data 
# because this allows for more efficient processing of large-scale data sets.
sparse.size <- object.size(pbmc.data)
sparse.size

# the ratio of the size of the dense matrix to the size of the sparse matrix
dense.size/sparse.size

### 2. QC and selecting cells for further analysis
# Seurat allows you to explore QC metrics and filter cells based on any user-defined criteria.

# Calculate and add the percentage information of mitochondrial genes to the Seurat object
# Low quality or dying cells often exhibit extensive mitochondrial contamination
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
# Including the number of RNA features (number of genes) per cell, 
#           the total count (expression),
#           and the percentage of mitochondrial genes.
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filter: the threshold is based on the first Vlnplot's results
# pbmc <- subset(pbmc,subset = nFeature_RNA <=2000 & nCount_RNA <= 5000 & percent.mt <= 5)
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter cells that have unique feature counts over 2500 
#                   or less than 200
#                   or have >5% mitochondrial counts
# on the basis of vlnplot's results
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### 3. Normalizing the data
# After removing unwanted cells from the dataset, the next step is to normalize the data
# Significance: Remove the impact of sequencing depth
# Principle: The count number of each gene in each cell is divided by the total count number of the cell
#            then multiplied by the factor (10000)
#            then log(n+1) transformed
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

### 4. Identification of highly variable features (feature selection) 识别高变基因
# return 2000 features per dataset that will be used in downstream analysis
# selection.method: vst, mvp, disp
# vst: 基于方差稳定性
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Selecting genes that are highly variable helps capture variability 
# between samples and thus better distinguish cell types or states.

# Identify the 10 most highly variable genes 
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### 5. Scaling the data
# scale all the genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

### 6. Perform linear dimensional reduction
# Only consider the previously identified highly variable genes 
# (obtained through the VariableFeatures function) for PCA dimensionality reduction.
# PCA reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

# DimHeatmap allows for easy exploration of the primary sources of heterogeneity
# and can be useful when trying to decide which PCs to include to further downstream analyses
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)



### 7. Determine the ‘dimensionality’ of the dataset
# Elbow
# In this example, we might have been justified in choosing anything between PC 7-12 as a cutoff.
ElbowPlot(pbmc)
# x-axis represents the number of the principal components
# y-axis represents the variance explanation rate

### 8. Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
# > head(Idents(pbmc), 5)
# AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 AAACCGTGTATGCG-1 
# 2                3                2                1                6 
# Levels: 0 1 2 3 4 5 6 7 8
# There are 9 clusters in this example

### 9. Run non-linear dimensional reduction (UMAP/t-SNE)
# note: avoid drawing biological conclusions solely on the basis of visualization techniques.

# UMAP based on the cluster results
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function 
# to help label individual clusters
DimPlot(pbmc, reduction = "umap")

### Save R data
saveRDS(pbmc, file = "pbmc_tutorial.rds")

### 10. Find differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
# FindMarkers: it identifies positive and negative markers of a single cluster (specified in ident.1)
#              compared to all other cells. If the ident.2 parameter is omitted or set to NULL
#              the FindMarkers function performs differential expression analysis between 
#              the specified ident.1 group and all other groups.
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
  
# find markers for every cluster compared to all remaining cells, 
# only.pos = TRUE: report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# avg_log2FC: Log value of the mean fold difference in expression between the two groups. 
#             Positive values indicate that the gene is more highly expressed in the first group.

# Seurat has several tests for differential expression which can be set with the "test.use" parameter 
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features.
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

### 11. Assigning the cell type identity to clusters
#
# Cluster ID    Markers	        Cell Type
# 0	            IL7R, CCR7	    Naive CD4+ T
# 1	            CD14, LYZ	      CD14+ Mono
# 2	            IL7R, S100A4	  Memory CD4+
# 3	            MS4A1	          B
# 4	            CD8A	          CD8+ T
# 5	            FCGR3A, MS4A7	  FCGR3A+ Mono
# 6	            GNLY, NKG7	    NK
# 7	            FCER1A, CST3	  DC
# 8	            PPBP	          Platelet

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", 
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot


###########################
# Mar. 20, 2024

### Analysis of Sequence-based Spatial Data

### libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ape)

### Data
# Mouse brain spatial expression split among across datasets.
InstallData("stxBrain")

brain <- LoadData("stxBrain", type = "anterior1")

### Data pre-processing

# check
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# SCTransform 构建基因表达的正则化负二项式模型
# 对数据进行标准化，检测高方差特征，并存储在SCT assay中
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

### Gene expression visualization
# Gene Hpca is a strong hippocampus marker
# Gene Ttr is a marker of the choroid plexus
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

plot <- SpatialFeaturePlot(brain, features = c("Ttr")) + 
  theme(legend.text = element_text(size = 0),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(1, "cm")
)
plot

# pt.size.factor: 点的大小 default = 1.6
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)

# alpha: 透明度 default = c(1,1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

### Dimensionality reduction, clustering, and visualization
# the same workflow as the scRNA-seq analysis
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)

# overlaid on the image
p2 <- SpatialDimPlot(object = brain, label = TRUE, label.size = 3)
p1 + p2

# Use cell.highlight to demarcate particular cells of interest on a SpatialDimPlot()
SpatialDimPlot(
  brain, 
  cells.highlight = CellsByIdentities(
    object = brain, 
    idents = c(2, 1, 4, 3, 5, 8)
  ), 
  facet.highlight = TRUE, ncol = 3
)

### Interactive plotting
# interactive parameter for Dimplot and FeaturePlot
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

# LinkedDimPlot links the UMAP representaion to the tissue image
# also allows interactive selection
LinkedDimPlot(brain)

### Identification of Spatially Variable Features

# Perform differential expression based on pre-annotated anatomical regions within the tissue
# which may be determined either from unsupervised clustering or prior knowledge.
# Find markers (DEG) for identity classes
# ident.1: 用于定义marker的identity class
# ident.2: 用于比较的第二个identity class; if NULL, use all other cells for comparison
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
# pct.1: The percentage of cells where the gene is detected in the first group
# pct.2: The percentage of cells where the gene is detected in the second group

SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:2],
                   alpha = c(0.1, 1), ncol = 1)

# Alternative: search for features exhibiting spatial patterning in the absence of pre-annotation
brain <- FindSpatiallyVariableFeatures(
  brain, assay = "SCT",
  features = VariableFeatures(brain)[1:1000],
  selection.method = "moransi"
)
# Moran's I: 莫兰指数 range from -1.0 to 1.0
# >0:正相关; <0:负相关; =0:随机; 绝对值越大，差异越大
# P: Probability, P很小表示所观察的空间不太可能产生随机过程
# Z: 标准差倍数，反映一个数据集的离散程度

#top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
#SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))

cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

# crop: crop the plot in to focus on points plotted. Set to FALSE to show entire image
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

### Intergration with single-cell data
allen_reference <- readRDS("allen_cortex.rds")

# Normalize
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# After subsetting we re-normalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# the annotation is stored in the "subclass" column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

cortex <- FindSpatiallyVariableFeatures(
  cortex, assay = "predictions",
  selection.method = "moransi",
  features = rownames(cortex),
  r.metric = 5,
  slot = "data"
)

#####################################
# Mar. 

### Analysis of Image-based Spatial Data

### libraries
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)


