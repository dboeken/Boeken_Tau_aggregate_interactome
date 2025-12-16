import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp

from loguru import logger

logger.info('Import OK')

input_path = 'results/preprocessed/normalised_summary.csv'
output_folder = 'results/preprocessed/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
 
# Read in dataset
clean_data = pd.read_csv(input_path)
clean_data.drop([col for col in clean_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)


# =============Calculate abundance ratio Disease / Control=============
control_average = clean_data[clean_data['disease'] == 'CRL'].copy()
control_average = control_average.groupby(['pulldown', 'antibody', 'pulldown-type', 'Protein IDs', 'Gene names', 'Protein names', 'Number of proteins']).mean()[['ref-corr_abundance', 'tau-corr_abundance']].reset_index()
control_average.rename(columns={'ref-corr_abundance': 'ref-corr_abundance-CRL', 'tau-corr_abundance': 'tau-corr_abundance-CRL'}, inplace=True)

# Merge control averages back onto original df
ratios = pd.merge(clean_data, control_average, on=['pulldown', 'antibody', 'pulldown-type', 'Protein IDs', 'Gene names', 'Protein names', 'Number of proteins'], how='outer')
ratios['RC_abundance-ratio'] = ratios['ref-corr_abundance'] / ratios['ref-corr_abundance-CRL']
ratios['log_RC_abundance-ratio'] = np.log2(ratios['RC_abundance-ratio'] )
ratios['TC_abundance-ratio'] = ratios['tau-corr_abundance'] / ratios['tau-corr_abundance-CRL']
ratios['log_TC_abundance-ratio'] = np.log2(ratios['TC_abundance-ratio'] )


# Calculate CTRL vs AD one-sample t-test
stats = []
for group, df in ratios[ratios['disease'] == 'AD'].groupby(['Protein IDs', 'Gene names', 'Protein names', 'Number of proteins', 'pulldown', 'antibody', 'pulldown-type']):
    df
    if len(df['log_RC_abundance-ratio'].dropna()) < 3:
        continue
    
    tval_1, pval_1 = ttest_1samp(
        df['log_RC_abundance-ratio'].dropna().tolist(),
        popmean=0
    )
    tval_2, pval_2 = ttest_1samp(
        df['RC_abundance-ratio'].dropna().tolist(),
        popmean=1
    )
    stats.append([*group, pval_1, pval_2, df['log_RC_abundance-ratio'].dropna().mean(), df['RC_abundance-ratio'].dropna().mean()])
summary = pd.DataFrame(stats, columns=['Protein IDs', 'Gene names', 'Protein names', 'Number of proteins', 'pulldown', 'antibody', 'pulldown-type', 'logRC-AR_pval', 'RC-AR_pval', 'log_RC_abundance-ratio', 'RC_abundance-ratio'])

summary['-log10_logRC-AR_pval'] = -np.log10(summary['logRC-AR_pval'])
summary['-log10_RC-AR_pval'] = -np.log10(summary['RC-AR_pval'])


# Add categories based on fold-change and significance 
summary['category'] = ['sig' if ((abs(fc) > 1) and (sig > 1.3)) else 'ns' for fc, sig in summary[['log_RC_abundance-ratio', '-log10_logRC-AR_pval']].values ]

summary['category'].value_counts()


# ====================Save to CSV====================
summary.to_csv(f'{output_folder}ratio_summary.csv')

# Collect significant proteins
sigs = summary[summary['category'] == 'sig'].copy()
# Collect only last entry for the protein groups --> need to come back and decide if this is the strategy we want to stick with
sigs['proteins'] = sigs['Protein IDs'].str.split(';').str[-1]

sig_list = {}
for (antibody, pulldown), df in sigs.groupby(['antibody', 'pulldown-type']):
    sig_list[f'{antibody}_{pulldown}_all'] = df['proteins'].tolist()
    sig_list[f'{antibody}_{pulldown}_up'] = df[df['log_RC_abundance-ratio'] > 0]['proteins'].tolist()
    sig_list[f'{antibody}_{pulldown}_down'] = df[df['log_RC_abundance-ratio'] < 0]['proteins'].tolist()

sig_list = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in sig_list.items() ]))

sig_list.to_csv(f'{output_folder}significant_summary.csv')

#  Add PERCEPT scaling

def percept(m0, m1, F, p):
    return m0 + ((m0 - m1) / -(F**p))

scaled = []
for (protein, pulldown), df in summary.groupby(['Protein IDs', 'pulldown']):
    
    df['percept_RC-AR'] = percept(m0=0, m1=df['log_RC_abundance-ratio'].values, F=10*4, p=df['logRC-AR_pval'].values)
    
    scaled.append(df)

scaled = pd.concat(scaled)
scaled.to_csv(f'{output_folder}percept_summary.csv')
