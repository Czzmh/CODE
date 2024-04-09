#####################################
# Mar. 25, 2024

### libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(SingleR)
library(celldex)
library(pheatmap)

# Data
data <- LoadH5Seurat("Data/SS200000116BR_E6.bin200.h5seurat")
# data is already Seurat object
head(data@meta.data, 5)
#                orig.ident    x     y
# 13743895370950     sample 3200 23750
# 20401094664900     sample 4750  8900
# 22333829957250     sample 5200 18050
# 31138512919850     sample 7250 23850
# 32212254743550     sample 7500 23550

# Calculate nCount_RNA and nFeature_RNA
nCount_RNA <-  colSums(GetAssayData(data)) # 计算每个细胞的RNA读数计数
nFeature_RNA <-  colSums(GetAssayData(data) > 0) # 计算每个细胞检测到的唯一基因数
data@meta.data$nCount_RNA <- nCount_RNA
data@meta.data$nFeature_RNA <- nFeature_RNA
head(data@meta.data, 5)

# MT
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
head(data@meta.data, 5) 

saveRDS(data, "RData/E6_bin200.rds")

data <- readRDS("RData/E6_bin200.rds")

# VlnPlot (save as VlnPlot1)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filter (based on the previous Vlnplot)
data <- subset(data, subset = nFeature_RNA >= 2500 & nFeature_RNA <= 10000 & percent.mt <= 10 & percent.mt >= 1)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

# 识别高变基因
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes 
top10 <- head(VariableFeatures(data), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scale
all.genes <- rownames(data)
data <- ScaleData(data)

### 6. Perform linear dimensional reduction
data <- RunPCA(data, features = VariableFeatures(object = data))

# Examine and visualize PCA results a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca") + NoLegend()

DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow
# In this example, we might have been justified in choosing anything between PC 7-10 as a cutoff.
ElbowPlot(data)


data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
head(Idents(data), 5)

# UMAP based on the cluster results
data <- RunUMAP(data, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function 
# to help label individual clusters
DimPlot(data, reduction = "umap")

VlnPlot(data, features = c("AL391095.2", "IGHM"))

### Save R data
saveRDS(data, file = "seurat_result.rds")

data <- readRDS("seurat_result.rds")

data.markers <- FindAllMarkers(data, only.pos = TRUE)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

### Annotation
### SingleR

load("Data/HumanPrimaryCellAtlas_hpca.se_human.RData")
load("Data/BlueprintEncode_bpe.se_human.RData")

data <- readRDS("seurat_result.rds")
data

meta <- data@meta.data
head(meta)
# UMAP
DimPlot(data, reduction = "umap", label = TRUE)

data_for_SingleR <- GetAssayData(data, slot="data")
data.hesc <- SingleR(test = data_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
data.hesc

# table of seurat and singleR
table(data.hesc$labels,meta$seurat_clusters)
