import pandas as pd
import NaiveDE
import SpatialDE
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Set data path
data_path = "D:/Desktop/MGI/CODE/Data/"

data = sc.read_h5ad(data_path + "SS200000116BR_E6.bin200.h5ad")

data_df = data.to_df()
colsum = data_df.sum(axis = 1)
data_df['total_counts'] = colsum

# pre-processing
sc.pp.filter_cells(data, min_genes = 500)  # filter低于200个genes的细胞
sc.pp.filter_cells(data, max_genes = 10000)
sc.pp.filter_genes(data, min_cells = 3) # filter在少于3个细胞中表达的genes
print(data.shape) # shape

# normalize
# target_sum: 每个样本的目标总表达量，归一化为1e4
sc.pp.normalize_total(data, target_sum=1e4)

sc.pp.log1p(data)
var_genes = sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.25, inplace = False)

# Get highly variable genes
hv_genes = var_genes[var_genes["highly_variable"] == True].index
print(hv_genes)

# Get the index of highly variable genes
hv_genes_index = np.where(var_genes["highly_variable"])[0]
print("\n", hv_genes_index)

data = data[data.obs_names, data.var_names[hv_genes_index]]

counts = data.to_df()
coords = data.obs[['x', 'y', 'total_counts']]

norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(coords, norm_expr.T, 'np.log(total_counts)').T

X = coords[['x','y']].to_numpy()
sample_resid_expr = resid_expr.sample(n=1000, axis=1, random_state=1)

results = SpatialDE.run(X, sample_resid_expr)