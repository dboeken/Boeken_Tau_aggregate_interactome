import os
import pandas as pd
import numpy as np

from src.utilities.statistical_tests import apply_enrichment

from loguru import logger
logger.info('Import OK')

input_path = 'results/preprocessed/RC_significant_summary.csv'
background_path = 'results/preprocessed/identified_background.csv'
output_folder = 'results/RC_go_enrichment/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# -----------------Read in standard components-----------------
# cleaned ratio summary
proteins = pd.read_csv(f'{input_path}')
proteins.drop([col for col in proteins.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# raw data as background
background = pd.read_csv(f'{background_path}')
background.drop([col for col in background.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
background = background['0'].tolist()


# ----------------Perform Panther enrichment test---------------- 

# perform enrichment test
enrichment = apply_enrichment(proteins, searches=None, obo_path='resources/bioinformatics_databases/PANTHERGOslim.obo', organism='9606', refOrganism='9606', enrichmentTestType='FISHER', correction='FDR', min_proteins=2, reference=background)

# Save all to excel
enrichment.to_csv(f'{output_folder}enrichment_complete.csv')
