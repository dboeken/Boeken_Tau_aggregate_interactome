import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import zscore
import numpy as np


input_path = 'results/preprocessed/normalised_summary.csv'
output_folder = 'results/plot_heatmap/'


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

import matplotlib
font = {'family': 'normal',
        'weight': 'normal',
        'size': 12}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 300

# Read in data
summary = pd.read_csv(input_path)
summary.drop([col for col in summary.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Exclude the sample "REF-0"
summary = summary[summary['sample'] != 'REF-0']

experiments = [col for col in summary['pulldown'].unique() if 'Beads' in col]

exp_data = summary[summary['pulldown'].isin(experiments)].copy()
exp_data['key'] = exp_data['pulldown'] + '_' + exp_data['sample']

# Pivot the data to have proteins as rows and samples (patients) as columns
heatmap_data = exp_data.pivot_table(index='Protein IDs', columns='key', values='ref-corr_abundance')

# Sort the columns by disease status (group patients by disease)
antibodies = ['AT8', 'T181', 'HT7', 'TOMA', 'MC1']
disease_states = ['AD', 'CRL']
donors = ['1', '2', '3', '4', '5']

sample_order = [f'{pulldown}-Beads_{disease}-{donor}' for pulldown in antibodies for disease in disease_states for donor in donors if f'{pulldown}-Beads_{disease}-{donor}' in exp_data['key'].unique()]

heatmap_data = heatmap_data[sample_order].copy()


heatmap_data = heatmap_data.dropna(how='any')

# Standardise (Z-score normalisation) across proteins 
norm_heatmap_data = heatmap_data.copy()
norm_heatmap_data[norm_heatmap_data.columns.tolist()] = zscore(heatmap_data, axis=1)




norm_heatmap_data = norm_heatmap_data.dropna()

# Perform hierarchical clustering
linkage_matrix = linkage(norm_heatmap_data, method='ward')

# Prune the tree by setting a threshold for the clustering distance
max_d = 10  
clusters = fcluster(linkage_matrix, max_d, criterion='distance')

sns.clustermap(
    norm_heatmap_data,
    row_cluster=True,  # Cluster proteins
    col_cluster=False,  # Cluster patients
    cmap='magma',  
    # standard_scale=0,  # Normalise each row independently
    figsize=(10, 26),  
    dendrogram_ratio=(0.2, 0.1),  #
    xticklabels=True,
    yticklabels=False
)

plt.savefig(f'{output_folder}heatmap2.svg')
plt.show()

