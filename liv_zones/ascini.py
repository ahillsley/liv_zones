import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from cellpose import utils
from cellpose.io import imread


def load_and_filter_cells(path):
    all_props = pd.read_csv(f'{path}/average_properties_per_cell.csv')
    
    # filter out bad cells based on thresholds
    
    
    return filtered_props
    

def plot_properties(
        path,
        colors=['steelblue', 'coral'],
        ):
    # make a new folder and save all individual plots within
    
    try:
        os.mkdir(f'{path}/ascini_trends/')
    except:
        print('file already exists, will rewrite over old plots')
    
    skip_props =[
        'Unnamed: 0',
        'area',
        'centroid-0',
        'centroid-1',
        'label',
        'ascini_position'
        ]
    
    all_props = pd.read_csv(
        f'{path}/average_properties_per_cell.csv').sort_values('ascini_position')
    
    position = all_props['ascini_position']
    
    for prop in all_props.keys():
        if prop in skip_props:
            continue
        
        single_prop = all_props[prop]
        
        trend = single_prop.rolling(20, min_periods=10, center=True).mean()
        
        plt.scatter(position, single_prop, c=colors[0])
        plt.plot(position, trend, c=colors[1], linewidth=3)
        plt.xlabel('Ascinus Position')
        plt.ylabel(prop)
        plt.savefig(
            f"{path}/ascini_trends/{prop}.png", dpi=300, bbox_inches="tight"
        )
        plt.figure()
        
        
    


def plot_cell(
    image_path,
    mask_path,
    cell_num,
    ):

    image = imread(image_path)
    cell_mask = np.load(f"{mask_path}/cell_mask.npy")

    bbox = _isolate_cell(cell_mask, cell_num)

    cell_mask_roi = cell_mask[bbox[0] : bbox[1], bbox[2] : bbox[3]] == cell_num

    mito_mask = np.load(f"{mask_path}/mito_mask.npy")
    mito_outline = utils.masks_to_outlines(mito_mask)
    mito_mask[mito_outline == 1] = 0

    lipid_mask = np.load(f"{mask_path}/lipid_droplet_mask.npy")
    lipid_outline = utils.masks_to_outlines(lipid_mask)
    lipid_mask[lipid_outline == 1] = 0

    actin_raw_roi = image[0, bbox[0] : bbox[1], bbox[2] : bbox[3]]
    mito_raw_roi = image[1, bbox[0] : bbox[1], bbox[2] : bbox[3]] * cell_mask_roi
    lipid_raw_roi = image[2, bbox[0] : bbox[1], bbox[2] : bbox[3]] * cell_mask_roi

    mito_mask_roi = mito_mask[bbox[0] : bbox[1], bbox[2] : bbox[3]] * cell_mask_roi
    lipid_mask_roi = lipid_mask[bbox[0] : bbox[1], bbox[2] : bbox[3]] * cell_mask_roi

    fig, axs = plt.subplots(2, 3)
    axs[0, 0].imshow(actin_raw_roi)
    axs[0, 1].imshow(mito_raw_roi)
    axs[0, 2].imshow(lipid_raw_roi)

    axs[1, 1].imshow(mito_mask_roi)
    axs[1, 2].imshow(lipid_mask_roi)
    [axi.set_axis_off() for axi in axs.ravel()]

    plt.savefig(
        f"{mask_path}/cell_{cell_num}_segmentation.png", dpi=300, bbox_inches="tight"
    )
    return


def plot_ascinus_annotated(path):
    # plot an ascinus with each cell labeled with its cell id
    cell_mask = np.load(f"{path}/cell_mask.npy")
    outline = utils.masks_to_outlines(cell_mask)
    cell_mask[outline == 1] = 0
    cell_mask[cell_mask != 0] += 50

    cell_props = pd.read_csv(f"{path}/average_properties_per_cell.csv")
    plt.imshow(cell_mask)
    plt.axis("off")
    for cell in cell_props["label"]:
        plt.text(
            x=cell_props["centroid-1"][cell - 1],
            y=cell_props["centroid-0"][cell - 1],
            s=f"{cell}",
            size=4,
            horizontalalignment="center",
            fontweight="semibold",
        )

    plt.savefig(f"{path}/labled_cells.png", dpi=300, bbox_inches="tight")
    return


def _isolate_cell(cell_mask, cell_num):
    a = np.where(cell_mask == cell_num)
    return [np.min(a[0]), np.max(a[0]), np.min(a[1]), np.max(a[1])]
