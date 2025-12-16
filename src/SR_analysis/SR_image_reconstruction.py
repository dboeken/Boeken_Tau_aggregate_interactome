# ## Manually export images using the picasso functions to render oversampled images with optional precision based blur.
from configparser import Interpolation
import os
import numpy as np
import pandas as pd
from picasso import render, io
import matplotlib.pyplot as plt
from skimage.io import imread
import matplotlib.patches as mpatches
import napari
import numpy as np
from skimage import data, transform

input_parameters = 'results/super-res/initial_cleanup/slide_parameters.csv'
input_folder = 'results/super-res/localisation/'
output_folder = 'results/super-res/rendering/'
oversampling = 8

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Compile list of files to be processed
parameters = pd.read_csv(input_parameters)
parameters.drop([col for col in parameters.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
parameters.dropna(subset=['keep'], inplace=True)

# Read in localisations
for filepath, well_info in parameters[['file_path', 'well_info']].values:
    filepath
    try:
        locs, info = io.load_locs(f'{input_folder}{well_info}_locs_corr_filt.hdf5')
        original = imread(filepath)
    except:
        continue

    # Get minimum / maximum localizations to define the ROI to be rendered 
    x_min = np.min(locs.x)    
    x_max = np.max(locs.x)
    y_min = np.min(locs.y)
    y_max = np.max(locs.y)

    viewport =  (y_min, x_min), (y_max, x_max)
    oversampling = 8
    len_x, image = render.render(locs, viewport = viewport, oversampling=oversampling, blur_method='smooth')
    plt.imsave(f'{output_folder}{well_info}_smooth.png', image, cmap='hot', vmax=10)

    fig, axes = plt.subplots(1, 2, figsize=(30, 15))
    axes[0].imshow(np.max(original, axis=0), cmap='Greys_r')
    axes[0].set_title('Original')
    axes[1].imshow(image, cmap='hot', vmax=8)
    axes[1].set_title('Super-resolved')
    plt.savefig(f'{output_folder}{well_info}_sr.png')
    plt.show()



# Cutom ROI with higher oversampling for a single image:
# Open necessary components
well = ''
filepath, well_info = parameters[parameters['well_info'] == well][['file_path', 'well_info']].values[0]
locs, info = io.load_locs(f'{input_folder}{well_info}_locs_corr_filt.hdf5')
original = imread(filepath)
original_image = transform.rescale(np.max(original, axis=0), 8, preserve_range=True)

viewport =  (np.min(locs.y), np.min(locs.x)), (np.max(locs.y), np.max(locs.x))
viewport =  (0, 0), (512, 512)
oversampling = 8
len_x, sr_image = render.render(locs, viewport = viewport, oversampling=oversampling, blur_method='smooth')
# ---------


with napari.gui_qt():
    # add the image
    viewer = napari.Viewer()
    viewer.add_image(original_image, name='Original', contrast_limits=(2, 65555))
    viewer.add_image(sr_image, name='Super resolved')
    shapes = viewer.add_shapes()
napari.run()

#[ystart, ystop, xstart, xstop],
viewports = [[int(ystart/oversampling), int(ystop/oversampling), int(xstart/oversampling), int(xstop/oversampling)] for [[ystart, xstart], [_, _], [ystop, xstop], [_, _]] in shapes.data]

zoomsampling = 20
palette = {0: 'orange', 1: 'red', 2: 'rebeccapurple'}
for x, chunk in enumerate(viewports):
    ystart, ystop, xstart, xstop = chunk
    len_x, zoom_image = render.render(locs, viewport = ((ystart, xstart), (ystop, xstop)), oversampling=zoomsampling, blur_method='smooth')
    plt.imsave(f'{output_folder}{well_info}_zoom_{x}_smooth.png', zoom_image, cmap='hot', vmax=10)
    
    fig, axes = plt.subplots(1, 2, figsize=(30, 15))
    axes[0].imshow(np.max(original[:, ystart: ystop, xstart: xstop], axis=0), cmap='Greys_r', interpolation='none')
    axes[0].set_title('Original')
    axes[0].axis('off')
    axes[1].imshow(zoom_image, cmap='hot', vmax=10)
    axes[1].set_title('Super-resolved')
    axes[1].axis('off')
    plt.savefig(f'{output_folder}{well_info}_zoom_{x}_sr.png')
    plt.show()


