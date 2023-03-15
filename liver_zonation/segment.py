from cellpose import models
from cellpose.io import imread
from cellpose import utils
import numpy as np
from scipy.ndimage import distance_transform_bf


def process_mito(img_path, mito_model, save_path):
    image = imread(img_path)
    mito = image[1,:,:]
    mito_masks = mito_model.eval(mito, flow_threshold=0.4, cellprob_threshold=0)
    
    np.save(f'{save_path}mito_mask.npy', mito_masks[0])
    return mito_masks[0]

def process_cells(img_path, cell_model, save_path):
    image = imread(img_path)
    actin = image[0,:,:]
    cell_masks = cell_model.eval(actin, flow_threshold=0.4, cellprob_threshold=0)
    
    np.save(f'{save_path}cell_mask.npy', cell_masks[0])
    return cell_masks[0]

def process_lipid_droplets(img_path, lipid_model, save_path):
    image = imread(img_path)
    ld = image[2,:,:]
    ld_masks = lipid_model.eval(ld, flow_threshold=0.4, cellprob_threshold=0)
    
    np.save(f'{save_path}ld_mask.npy', ld_masks[0])
    return ld_masks[0]

def portal_distance(save_path):
    pv_image = imread(f'{save_path}pv_distance.tif')
    pv_distance = distance_transform_bf(pv_image)
    np.save(f'{save_path}pv_distance.npy', pv_distance)
    
    return pv_distance

def central_distance(save_path):
    cv_image = imread(f'{save_path}cv_distance.tif')
    cv_distance = distance_transform_bf(cv_image)
    np.save(f'{save_path}cv_distance.npy', cv_distance)
    
    return cv_distance

def cell_edge_distance(save_path, cell_mask=None):
    if cell_mask is None:
        cell_mask = np.load(f'{save_path}cell_mask.npy')
        
    dist_transform = utils.distance_to_boundary(cell_mask)
    
    np.save(f'{save_path}cell_boundry_distance.npy', dist_transform)
    return

def outlines(image):
    return utils.masks_to_outlines(image)
