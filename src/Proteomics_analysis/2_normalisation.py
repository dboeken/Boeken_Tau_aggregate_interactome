import os, re
import numpy as np
import pandas as pd

from loguru import logger

logger.info('Import OK')


def norm_control_plex(peptides, sample_cols, pooled_col, standard_vals):
    """ Normalises per peptide according to control channel. Returns df with updated sample cols and removes control channel"""
    pooled_factor = pd.merge(pd.DataFrame(peptides[pooled_col]), pd.DataFrame(standard_vals), left_index=True, right_index=True)
    pooled_factor['pooled_factor'] = pooled_factor.iloc[:, 0] / pooled_factor.iloc[:, 1]
    
    # normalise to control channel
    peptides[sample_cols] = peptides[sample_cols].multiply(
        pooled_factor['pooled_factor'], axis=0)

    return peptides


def noncys_ci_calculator(noncys_peptides, sample_cols):
    """Calculates relevant info per protein from noncys peptides, including overall mean, SD
    and the per-channel mean (used for cys/noncys calculation)"""

    # for each protein, collect all noncys values and determine mean + SD
    noncys_cis = {}
    for protein, df in noncys_peptides.groupby(['Proteins']):
        values = df[sample_cols].values.flatten()
        # Remove NaN values for calculations
        values = values[~np.isnan(values)]
        noncys_cis[protein] = {'df': df,
                               'noncys_means': df[sample_cols].mean(),
                               'num_peptides': df.shape[0],
                               'values': values,
                               'pop_mean': np.mean(values),
                               'pop_stdev': np.std(values),
                               'num_vals': len(values)}
    return noncys_cis


def med_normalise(peptides, control_plex):
    """Calculates correction factor from median of column 1, then normalises all other channels to this factor"""

    medians = peptides.median()
    control_factor = medians[control_plex] / medians

    # normalise to control channel
    peptides = peptides * control_factor

    return peptides


#---------------------------------------------------------------------------------------------

if __name__ == "__main__":

    input_folder = 'results/initial_cleanup/'
    output_folder = 'results/preprocessed/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_data = pd.read_csv(f'{input_folder}Compiled_Proteins.csv')
    raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
    
    # Collect background
    pd.DataFrame(raw_data['Protein IDs'].str.split(';').str[-1].unique()).to_csv(f'{output_folder}identified_background.csv')

    clean_data = raw_data.copy()
    
    # ===============Normalise to reference channel (131)===============
    # Normalise according to sum of the reference channel - 
    # should correspond to loading equivalent total protein
    ref_channels = clean_data[clean_data['disease'] == 'REF'].copy()
    # Determine scaling factor based on maximum total abundance in a given exp
    sum_refs = ref_channels.groupby('pulldown').sum()['intensity'].reset_index()
    sum_refs['corr_factor'] = sum_refs['intensity'].mean() / sum_refs['intensity'] 
    corr_factor = dict(sum_refs[['pulldown', 'corr_factor']].values)

    # Return scaling factor to original dataset
    clean_data['ref_corr-factor'] = clean_data['pulldown'].map(corr_factor)
    clean_data['ref-corr_abundance'] = clean_data['intensity'] * clean_data['ref_corr-factor']

    # ==================Normalise to Tau concentration==================
    tau_vals = clean_data.copy()
    # Determine scaling factor based on maximum total abundance in a given exp
    tau_vals['ref-corr_abundance'] = tau_vals['ref-corr_abundance'].astype(float)
    sum_refs = tau_vals[tau_vals['Protein IDs'] == 'P10636'].groupby(['pulldown', 'sample']).mean()['ref-corr_abundance'].reset_index()
    sum_refs['corr_factor'] = sum_refs['ref-corr_abundance'].mean() / sum_refs['ref-corr_abundance'] 
    sum_refs['key'] = [f'{exp}_{val}' for exp, val in sum_refs[['pulldown', 'sample']].values]
    corr_factor = dict(sum_refs[['key', 'corr_factor']].values)

    # Return scaling factor to original dataset
    clean_data['tau_corr-factor'] = [corr_factor[f'{exp}_{val}'] for exp, val in clean_data[['pulldown', 'sample']].values]
    clean_data['tau-corr_abundance'] = clean_data['ref-corr_abundance'] * clean_data['tau_corr-factor']



    # =========Remove AD-2 outlier=========
    
    clean_data = clean_data[clean_data['sample'] != 'AD-2'].copy()

    clean_data.to_csv(f'{output_folder}normalised_summary.csv')
    