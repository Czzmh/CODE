# Mar. 14, 2024

# scPred is a general method to classify cells based on a low-dimensional representation of gene expression (e.g. PCA)
# It uses machine learning methods, such as random forests and support vector machines, to classify and annotate single-cell RNA sequencing data.

### install packages
# BiocManager::install("scPred")
devtools::install_github("powellgenomicslab/scPred")

### libraries
library("scPred")
library("Seurat")
library("magrittr")
library("dplyr")

# Use the PBMCs from one individual to build cell classifiers for the populations of interest.
# Then, apply these models to an indepent dataset of PBMCs from another independent individual.

reference <- scPred::pbmc_1
query <- scPred::pbmc_2

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  # UMAP using the top 30 most variable PCs
  RunUMAP(dims = 1:30)

# The column cell_type contains the identity of each cell in the metadata slot. 
# Plot the UMAP and grouping the cells by cell type.
DimPlot(reference, group.by = 'cell_type', label = TRUE, repel = TRUE)

### Training classifiers with scPred
# Firstly, let’s get the feature space to train the classifiers. 
# By default, scPred will use all principal components. 

# getFeatureSpace will create a scPred object stored in the @misc slot. 
# This object will contained all required information to classify cells. 
reference <- getFeatureSpace(reference, "cell_type")

# Secondly, we train the classifiers for each cell using the trainModel function. 
# By default, scPred will use a SVM with a radial kernel.
reference <- trainModel(reference)

# Training probabilities for each cell in the reference data 
# can be accessed using the get_probabilities method
get_probabilities(reference) %>% head(5)

# Printing a scPred object will show for each cell type:
# number of cells
# number of features used to train the model
# prediction model
# performance metrics
get_scpred(reference)

# Plot with the probability distribution for each cell type
plot_probabilities(reference)
# According to the plot, we can observe an overall lower performance for cMono and ncMono
# Other models may show an better performance

### Cell classification
# An important requirement for classifying cells is using the same normalization method for both the reference and the query datasets.

# First, normalize the query dataset (cells to be classfied).
query <- NormalizeData(query)

query <- scPredict(query, reference)
#
# Error in `GetAssayData()`:
#   ! `assay` must be one of "RNA", not "data".
# Maybe because of the Seurat version









