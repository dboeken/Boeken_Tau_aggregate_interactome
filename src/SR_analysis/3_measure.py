import os
import pandas as pd
import numpy as np

from skimage import measure
from skan import Skeleton, summarize
from smma.utilities import measure_properties
# from smma.visualise import plot_branch_overview, plot_skeletonisation

from loguru import logger
logger.info('Import OK')

# root_path = open('raw_data/root_path.txt', 'r').readlines()[0]
root_path = ''


input_parameters = 'results/super-res/initial_cleanup/slide_parameters.csv'
input_folder = f'{root_path}results/super-res/localisation/'
output_folder = f'{root_path}results/super-res/measurements/'
spacing_nm = 107 / 8

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# read parameters
parameters = pd.read_csv(input_parameters)
parameters.drop([col for col in parameters.columns.tolist()
                if 'Unnamed: ' in col], axis=1, inplace=True)
parameters.dropna(subset=['keep'], inplace=True)

# If layouts not present
layouts = {
    '': '1',
    '': '2',
}
parameters['layout'] = parameters['file_path'].str.split(
    '/').str[-2].map(layouts)




# Measure clustered and smooth skeletonised localisations
properties = measure_properties(
    parameters=parameters, input_folder=input_folder, spacing_nm=spacing_nm, summarise=True, output_folder=None)
properties.to_csv(f'{output_folder}properties_compiled.csv')




