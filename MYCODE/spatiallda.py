import functools
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
import scipy
from sklearn.model_selection import train_test_split
import time
import tqdm
from anndata import AnnData
import scanpy as sc

# Spatial LDA imports
import spatial_lda.model

from featurization import featurize_tumors
from featurization import neighborhood_to_marker
from featurization import make_merged_difference_matrices
from visualization import plot_samples_in_a_row
from visualization import plot_one_tumor_all_topics
from visualization import plot_one_tumor_topic
from visualization import plot_topics_heatmap
from visualization import plot_adjacency_graph

# Parameters
TRAIN_SIZE_FRACTION = 0.9
N_PARALLEL_PROCESSES = 8
TRAIN_SIZE_FRACTION = 0.9 
N_TOPICS = 5 
DIFFERENCE_PENALTY = 250 
RETRAIN_MODEL = True

# Set data path
data_path = "D:/Desktop/MGI/CODE/"

# data
data = sc.read_h5ad(data_path + "Data/SS200000116BR_E6.bin200.h5ad")

# annotation data
anno_data = pd.read_csv(data_path + "Output/RCTD_results.csv", index_col=0)
anno_data.index = anno_data.index.rename('ID')
anno_data.reset_index(inplace=True)
anno_data.drop(['x', 'y'], axis=1, inplace=True)

# merge
data_obs = data.obs.copy()
data_obs.index = data_obs.index.rename('ID')
data_obs.reset_index(inplace=True)

data_obs['ID'] = data_obs['ID'].astype(str)
anno_data['ID'] = anno_data['ID'].astype(str)
merged_data = pd.merge(data_obs, anno_data, on='ID', how='inner')

# is_tumor
merged_data['is_tumor'] = np.where(merged_data['final_type'] == 'Tumor-Cholang', True, False)

# is_immune
immune_cell_types = {'B cell': True, 'T-NK': True, 'Macrophage': True}
merged_data['is_immune'] = np.where(merged_data['final_type'].isin(immune_cell_types.keys()), 
                                    True, False)

# tumor data
tumor_data = merged_data[(merged_data['is_tumor'] == True)|(merged_data['is_immune'] == True)]
tumor_data = tumor_data[(tumor_data['x'] > 5900) & (tumor_data['y'] > 5500) 
                        & (tumor_data['x'] < 25900) & (tumor_data['y'] < 25700)] # filter

# Dataframe to Anndata
tumor_data_ann = AnnData(obs=tumor_data, var=data.var)
aligned_X = data[data.obs.index.isin(tumor_data_ann.obs["ID"])].X
tumor_data_ann = AnnData(X=aligned_X.toarray(), obs=tumor_data, var=data.var)

features = pd.DataFrame(tumor_data_ann.X.astype(int), columns=tumor_data_ann.var.index)
features.index = [(1, idx) for idx in features.index]

patient = tumor_data_ann.obs
patient.reset_index(drop=True, inplace=True)

patient_dict = dict()
patient_dict[1] = patient


##############
with open(data_path + "Data/patient_dfs.pkl", 'rb') as f:
    patient_dfs = pickle.load(f)
for patient_id in patient_dfs.keys():
    df = patient_dfs[patient_id]
    df['combined_cluster_id'] = (df['immune_cluster'].fillna(0) + 
                                (df.cluster_id + 12).fillna(0))
    df.loc[df['combined_cluster_id'] == 0, 'combined_cluster_id'] = None
    df.loc[:, 'is_tumor'] = ~df['isimmune']
    patient_dfs[patient_id] = df

with open(data_path + "Data/tumor_marker_features.pkl", 'rb') as f:
    tumor_marker_features = pickle.load(f)
##############

_sets = train_test_split(features, test_size=1.-TRAIN_SIZE_FRACTION)
train_features, test_features = _sets
train_difference_matrices = make_merged_difference_matrices(train_features, patient_dict, 'x', 'y')

def make_plot_fn(difference_matrices):  
    def plot_fn(ax, tumor_idx, features_df, patient_dfs):
        plot_adjacency_graph(ax, tumor_idx, features_df, patient_dfs, difference_matrices)
    return plot_fn
_plot_fn = make_plot_fn(train_difference_matrices)

plot_samples_in_a_row(train_features, _plot_fn, patient_dict) #, tumor_set=[3])

# # train lda model
# spatial_lda_model = spatial_lda.model.train(train_features, 
#                                             train_difference_matrices, 
#                                             n_topics=N_TOPICS, 
#                                             difference_penalty=DIFFERENCE_PENALTY, 
#                                             verbosity=1,
#                                             n_parallel_processes=N_PARALLEL_PROCESSES,
#                                             n_iters=3,
#                                             admm_rho=0.1,
#                                             primal_dual_mu=2)

# # save model
# with open(data_path + "Model/spatial_lda_tmp.pkl", 'wb') as f:
#     pickle.dump(spatial_lda_model, f)

# load model
with open(data_path + "Model/spatial_lda_tumor_immune.pkl", 'rb') as f:
    complete_lda = pickle.load(f)

plot_samples_in_a_row(complete_lda.topic_weights, plot_one_tumor_all_topics, patient_dict)

for t in range(N_TOPICS):
  plot_samples_in_a_row(complete_lda.topic_weights.iloc[:, t], 
                        plot_one_tumor_topic, patient_dict, tumor_set=None)