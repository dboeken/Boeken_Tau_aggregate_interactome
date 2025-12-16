import os
import shutil
import pandas as pd
import numpy as np
import smma.src.utilities as utilities
import napari
import functools
from skimage import io

from loguru import logger

logger.info('Import OK')

input_folder = utilities.locate_raw_drive_files(
    input_path='raw_data/raw_data.txt')
output_folder = 'results/spot_detection/initial_cleanup/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ---------------Initial_cleanup---------------

image_ids = {
    '1': 'file name',    
}

compiled_parameters = []
for layout, image_id in image_ids.items():

    # Generate layout file
    utilities.slide_layout_to_map(
        image_id,
        slide_layout=f'raw_data/slide_layout_{layout}.xlsx',
        output_folder='raw_data/'
    )

    # Collect list of images
    file_list = utilities.find_files(
        input_folder=f'{input_folder}{image_id}/', file_type='.tif')


    # Prepare slide parameters
    parameters = utilities.slide_map_to_params(
        file_list,
        input_layout=f'raw_data/layout_{image_id}.csv'
    )

    # ------------Image quality check------------
    # Preview images to check dimensions, useability of individual images - generate a list of FOVs to ignore
    if not os.path.exists(f'{output_folder}quality_check_{layout}.csv'):
        quality_check = utilities.quality_check_tiled(parameters, frames =25)
        quality_check = quality_check[['well_info', 'keep']].copy()
        quality_check.to_csv(f'{output_folder}quality_check_{layout}.csv')
    else:
        quality_check = pd.read_csv(f'{output_folder}quality_check_{layout}.csv')

    parameters['layout'] = layout
    parameters['keep'] = parameters['well_info'].map(
        dict(quality_check[['well_info', 'keep']].values))

    compiled_parameters.append(parameters)

# save parameters to csv
parameters = pd.concat(compiled_parameters)
parameters.to_csv(f'{output_folder}slide_parameters.csv')
parameters.fillna(0).groupby('keep').count()
