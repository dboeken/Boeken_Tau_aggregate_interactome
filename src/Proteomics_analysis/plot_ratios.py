import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore



from loguru import logger

logger.info('Import OK')

input_path = 'results/preprocessed/ratio_summary.csv'
output_folder = 'results/plot_ratios/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
import matplotlib
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 16 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Read in data
summary = pd.read_csv(input_path)
summary.drop([col for col in summary.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)


palette = {
    'sig': '#F06449',
    'ns': '#687578',
}

# ===============Plot volcano plots===============
for (antibody, pulldown_type), df in summary.groupby(['antibody', 'pulldown-type']):
    df
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.scatterplot(
        data=df,
        x='log_RC_abundance-ratio',
        y='-log10_logRC-AR_pval',
        hue='category',
        palette=palette
    )
    plt.title(f'{antibody} {pulldown_type}')
    ax.axhline(1.3, linestyle='--', color='black')
    ax.axvline(1, linestyle='--', color='black')
    ax.axvline(-1, linestyle='--', color='black')
    ax.set(xlabel='Log$_2$ Abundance ratio', ylabel='- Log$_{10}$ $p$-value')
    plt.tight_layout()
    ax.set(xlim=(-4, 4))
    # plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.legend('')
    plt.savefig(f'{output_folder}{antibody}-{pulldown_type}_volcano.svg')
    plt.show()

