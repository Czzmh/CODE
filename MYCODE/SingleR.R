# Mar. 14, 2024

### install packages
# BiocManager::install("SingleR")
# BiocManager::install("celldex")

### libraries
library(SingleR)
library(celldex)
library(Seurat)
library(pheatmap)

### Download annotation database
# hpca.se <- HumanPrimaryCellAtlasData()
# hpca.se
# save(hpca.se, file = "HumanPrimaryCellAtlas_hpca.se_human.RData")
# 
# bpe.se <- BlueprintEncodeData()
# bpe.se
# save(bpe.se, file = "BlueprintEncode_bpe.se_human.RData")

load("Data/HumanPrimaryCellAtlas_hpca.se_human.RData")
load("Data/BlueprintEncode_bpe.se_human.RData")

# Load PBMC RDS in Seurat tutorial
pbmc <- readRDS('RData/pbmc_tutorial.rds')
pbmc

meta <- pbmc@meta.data
head(meta)

# UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE)

# SingleR annotation by built-in dataset
pbmc_for_SingleR <- GetAssayData(pbmc, slot="data")
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
pbmc.hesc

# table of seurat and singleR
table(pbmc.hesc$labels,meta$seurat_clusters)

# plot UMAP
pbmc@meta.data$labels <- pbmc.hesc$labels
DimPlot(pbmc, group.by = c('seurat_clusters', 'labels'), reduction = 'umap')

# Annotate using BP and HPCA database
pbmc_1 <- pbmc
pbmc_1.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP = bpe.se, HPCA = hpca.se),
                       labels = list(bpe.se$label.main, hpca.se$label.main))
table(pbmc_1.hesc$labels, meta$seurat_clusters)

pbmc_1@meta.data$lables <- pbmc_1.hesc$labels
DimPlot(pbmc_1, group.by = c('seurat_clusters', 'labels'), reduction = 'umap')

### Annotation diagnosis
plotScoreHeatmap(pbmc.hesc)

# Plot the distribution of deltas 
# (i.e., the gap between the assignment score for the assigned label 
# and those of the remaining labels) across cells assigned to each reference label.
plotDeltaDistribution(pbmc.hesc, ncol = 3)

# compare with cluster results
tab <- table(label = pbmc.hesc$labels,
             cluster = meta$seurat_clusters)
pheatmap(log10(tab + 10))

###
# pruneScores() function can be used to remove the low-quality annotation
# based on the delta value
