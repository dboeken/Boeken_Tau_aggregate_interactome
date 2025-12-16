import os
import pandas as pd
from loguru import logger

from src.utilities.databases import network_interactions, create_uniprot_xref

logger.info('Import OK')

input_path = 'results/preprocessed/ratio_summary.csv'
output_folder = 'results/protein_interactions/'
resource_folder = 'resources/bioinformatics_databases/'

confidence_threshold = 0.5

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# -----Read in peptides data for all proteins identified-----

raw_data = pd.read_csv(f'{input_path}')
raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
proteins = raw_data['Protein IDs'].unique().tolist()

# Map genes to STRING ids
string_map = create_uniprot_xref(resource_folder, tax_id='9606', gene_ids=[], id_type='STRING')
string_map = dict(string_map[['UniProtKB-AC', 'ID']].values)
string_ids = [string_map[gene] for gene in proteins if gene in list(string_map.keys()) ]

# Colect all network interactions - remember limit is 2000 genes
interactions = network_interactions(genes=string_ids, tax_id='9606', id_type='string')

# Add back original gene IDs based on mapping
inv_map = {v: k for k, v in string_map.items()}
interactions['Protein_A'] = interactions['stringId_A'].map(inv_map)
interactions['Protein_B'] = interactions['stringId_B'].map(inv_map)

interactions.to_csv(f'{output_folder}STRING_protein_interactions.csv')



