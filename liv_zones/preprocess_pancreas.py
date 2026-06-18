# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 09:52:31 2025

@author: adhikarir
"""

import numpy as np
from numpy import asarray
import cv2 as cv
from PIL import Image

import glob
import warnings
from scipy.ndimage import distance_transform_edt  # uninstall and re-install

warnings.filterwarnings("ignore")

from cellpose import models, utils, plot
from cellpose.io import imread


from liv_zones import organelle_model_pancreas as org_models
#from liv_zones.distance_to_veins import main as vein_dist
import os
import torch

# Always define save path not including the last /

def preprocessing(image_path, save_path, channels, feature_list=None):
    # function to run segmentation of all organelles and get all
    # distance transforms needed for post_processing

    if feature_list is None:
        feature_list = file_check(save_path)

    features = " ".join(feature_list)

    if "cell" in features:
        print("segmenting cells")
        cell_model = org_models.OrganelleModel("pancreas_cell")
        cell_mask = cell_model.segment(
            img_path=image_path,
            channel=channels["actin"],
            save=True,
            save_path=save_path,
        )

        corrected_cell_mask = remove_small_masks(save_path=save_path, save=True,)

    

    if "mito" in features:
        print("segmenting mitos")
        mito_model = org_models.OrganelleModel("pancreas_mito")
        mito_mask = mito_model.segment(
            img_path=image_path,
            channel=channels["mito"],
            save=True,
            save_path=save_path,

        )

    if "lipid" in features:
        print("segmenting lipid droplets")
        lipid_model_small = org_models.OrganelleModel("lipid_droplet_Small")
        lipid_mask_small = lipid_model_small.segment(

            img_path=image_path,
            channel=channels["lipid"],
            save=True,
            save_path=save_path,)
        
        lipid_model_medium = org_models.OrganelleModel("lipid_droplet_Medium")
        lipid_mask_medium = lipid_model_medium.segment(

            img_path=image_path,
            channel=channels["lipid"],
            save=True,
            save_path=save_path,)
        
        lipid_model_large = org_models.OrganelleModel("lipid_droplet_Large")
        lipid_mask_large = lipid_model_large.segment(
             img_path=image_path,
            channel=channels["lipid"],
            save=True,
            save_path=save_path,)
         
        #Using the Combining_lipid_droplets function combine the mask of 
        #small and large lipid droplets (Here also a .png 
        # of large, small and combine droplets masks is saved)
        
        lipid_mask_small_sm = f"{save_path}/lipid_droplet_Small_mask.npy"
        Small_droplets_array = np.load(lipid_mask_small_sm)
        #Small_droplets= cv.imread(lipid_mask_small_sm)
        #Small_droplets_array = asarray(Small_droplets)
        #SLD = Image.fromarray(Small_droplets_array)
        #SLD.save(f"{save_path}/small_lipid_droplets_mask.png")
        
        
        lipid_mask_large_lg = f"{save_path}/lipid_droplet_Large_mask.npy"
        Large_droplets_array_pre = np.load(lipid_mask_large_lg)
        Large_droplets_array = (Large_droplets_array_pre * (80)) 
        #Large_droplets= cv.imread(lipid_mask_large_lg)
        #Large_droplets_array = asarray(Large_droplets) #* 70 this multiplication helps visuslize the large droplets
        #LLD = Image.fromarray(Large_droplets_array) 
        #LLD.save(f"{save_path}/large_lipid_droplets_mask.png")
        
       
        lipid_mask_medium_md = f"{save_path}/lipid_droplet_Medium_mask.npy"
        Medium_droplets_array_pre = np.load(lipid_mask_medium_md)
        Medium_droplets_array = (Medium_droplets_array_pre * (10)) 

        #Generating a medium + small image by np.add
        Medium_and_Small_c = np.add(Small_droplets_array, Medium_droplets_array)
        np.save(f"{save_path}/lipid_droplet_mask_MS.npy", Medium_and_Small_c)

        Medium_and_Small_p = f"{save_path}/lipid_droplet_mask_MS.npy"
        Medium_and_Small = np.load(Medium_and_Small_p)

        #Generate the final mask containing small, medium ,and large lipid droplets
    
        lipid_mask = Combining_lipid_droplets(Medium_and_Small, Large_droplets_array,
            save=True,
            save_path=save_path, )
        #----------------------------------------------------------------------
     
    if "perox" in features:
        print("segmenting peroxisomes")
        peroxisome_model = org_models.OrganelleModel("pancreas_peroxisome")
        peroxisome_mask = peroxisome_model.segment(
            img_path=image_path,
            channel=channels["peroxi"],
            save=True,
            save_path=save_path,
        )
        
    if "large_perox" in features:
        print("segmenting large_peroxisomes")
        large_peroxisome_model = org_models.OrganelleModel("pancreas_large_peroxisome")
        large_peroxisome_mask = large_peroxisome_model.segment(
            img_path=image_path,
            channel=channels["peroxi"],
            save=True,
            save_path=save_path,
        )

    if "nuclei" in features:
        print("segmenting nuclei")
        nuclei_model = org_models.OrganelleModel("nuclei")
        nuclei_mask = nuclei_model.segment(
            img_path=image_path,
            channel=channels["nuclei"],
            save=True,
            save_path=save_path,
        )

    #if "central" in features or "portal" in features:
        #print("calculating distance to central and portal veins")
        #vein_distance = vein_dist(f'{save_path}/peroxisome_mask.npy', save_path)

    if "bound" in features:
        print("calculating distance to cell boundary")
        cell_edge_distance(save_path)

    return

# Funtion to add 2 masks from small and large lipid droplets

def Combining_lipid_droplets(image1,image2,save=True, save_path=""):

  
   All_lipid_droplets = np.add(image1,image2)
   
   if save is True:
       np.save(f"{save_path}/lipid_droplet_mask.npy", All_lipid_droplets)
    
       combine = Image.fromarray(All_lipid_droplets) 
       combine.save(f"{save_path}/lipid_droplets_mask.png")
       
   return All_lipid_droplets
    
  
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
def file_check(path):
    # check if any files are missing

    available_files = glob.glob(f"{path}/*")

    # Need to add mask files for any new organelles
    required_np_files = [
        "cell_mask.npy",
        "mito_mask.npy",
        "lipid_droplet_mask.npy",
        "peroxisome_mask.npy",
        "large_peroxisome_mask.npy",
        #"nuclei_mask.npy",
        #"central_dist.npy",
        #"portal_dist.npy",
        "boundary_dist.npy",
    ]

    missing = False
    to_process = []
    for file in required_np_files:
        file_path = f"{path}/" + file
        if file_path not in available_files:
            missing = True
            print(f"missing {file}")
            to_process.append(file)

    if missing is False:
        print("found all required files :)")
    print("-" * 30)

    return to_process

"""
def vein_distance(save_path, vein):
    pv_image = imread(f"{save_path}/{vein}v_distance.tif")
    pv_distance = distance_transform_edt(pv_image)
    np.save(f"{save_path}/{vein}v_distance.npy", pv_distance)

    return pv_distance
"""

def cell_edge_distance(save_path):
    cell_mask = np.load(f"{save_path}/cell_mask.npy")
    outlines = utils.masks_to_outlines(cell_mask)
    dist_transform = distance_transform_edt((outlines == False) * 1)

    np.save(f"{save_path}/boundary_dist.npy", dist_transform)
    return


def remove_small_masks(save_path, save=True):
   
    os.mkdir(f'{save_path}/small_cells_filtering/')

    img = np.load(f"{save_path}/cell_mask.npy")
    np.save(f"{save_path}/small_cells_filtering/pre_cell_mask.npy", img)
    To_PNG_1 = Image.fromarray(img) 
    To_PNG_1.save(f"{save_path}/small_cells_filtering/pre_cell_mask.png")
    torch.cuda.empty_cache()

    img2 = np.load(f"{save_path}/cell_mask.npy") 

    corrected = utils.fill_holes_and_remove_small_masks(img2, min_size= 1000) #min_size=110000)
    np.save(f"{save_path}/small_cells_filtering/corrected_cell_mask.npy", corrected)
   
    color_mask = plot.mask_rgb(corrected)
    color = Image.fromarray(color_mask) 
    color.save(f"{save_path}/small_cells_filtering/color_cell_mask.png")
   
    To_PNG_2 = Image.fromarray(corrected) 
    To_PNG_2.save(f"{save_path}/small_cells_filtering/corrected_cell_mask.png")
    torch.cuda.empty_cache()
    
    if save is True:
        cell_mask = np.load(f"{save_path}/small_cells_filtering/corrected_cell_mask.npy")
        np.save(f"{save_path}/cell_mask.npy", cell_mask)
        
        
    return 


if __name__ == "__main__":
    test_img_path = "../test_set/0-Actin_Dendra_LD-1_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif"
    path = "../test_set"
    a = file_check(path)
    preprocessing(test_img_path, path, a)
