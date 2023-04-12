import segment
import extract_properties as ex_prop
import visualize
from cellpose import models
from cellpose.io import imread
from cellpose import utils
import time


if __name__ == '__main__':
    
    save_paths = [
                    '../processed_imgs/0119_exp/ds1_m1_l2/0-image/',
                    '../processed_imgs/0119_exp/ds1_m1_l2/1-image/',
                    '../processed_imgs/0119_exp/ds1_m1_l4/0-image/',
                    '../processed_imgs/0119_exp/ds1_m1_l4/1-image/',
                    '../processed_imgs/0119_exp/ds1_m1_l4/2-image/',
                    '../processed_imgs/0119_exp/ds1_m3_l2/0-image/',
                    '../processed_imgs/0119_exp/ds1_m3_l3/0-image/',
                    '../processed_imgs/0119_exp/ds1_m3_l3/1-image/',
                    
                    '../processed_imgs/0119_exp/ds1_m3_l3/2-image/',
                    '../processed_imgs/0119_exp/ds2_m1_l1/0-image/',
                    '../processed_imgs/0119_exp/ds2_m1_l1/1-image/',
                    '../processed_imgs/0119_exp/ds2_m1_l1/2-image/',
                    '../processed_imgs/0119_exp/ds2_m1_l2/0-image/',
                    '../processed_imgs/0119_exp/ds2_m1_l2/1-image/',
                    '../processed_imgs/0119_exp/ds2_m1_l2/2-image/',
                    '../processed_imgs/0119_exp/ds2_m2_l2/0-image/',
                    '../processed_imgs/0119_exp/ds2_m2_l2/1-image/',
                    '../processed_imgs/0119_exp/ds2_m2_l2/2-image/',
                    '../processed_imgs/0119_exp/ds2_m3_l1/0-image/',
                    '../processed_imgs/0119_exp/ds2_m3_l1/1-image/',
                    '../processed_imgs/0119_exp/ds2_m3_l1/2-image/',
                    '../processed_imgs/0119_exp/ds3_m1_l1/2-image/',
                    '../processed_imgs/0119_exp/ds3_m1_l2/0-image/',
                    '../processed_imgs/0119_exp/ds3_m2_l1/0-image/',
                    '../processed_imgs/0119_exp/ds3_m2_l1/1-image/',   
                    '../processed_imgs/0119_exp/ds3_m2_l2/0-image/',
                    '../processed_imgs/0119_exp/ds3_m3_l1/0-image/',
                    '../processed_imgs/0119_exp/ds3_m3_l1/1-image/',
                    '../processed_imgs/0119_exp/ds3_m3_l1/2-image/',
                    '../processed_imgs/0119_exp/ds3_m3_l2/0-image/',
                    '../processed_imgs/0119_exp/ds3_m3_l2/1-image/',
                 ]
    
    images = [
              '../test-Dataset_1/m1/lobule2/0-Actin_Dendra_LD-1_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
               '../test-Dataset_1/m1/lobule2/1-Actin_Dendra_LD-2_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
               '../test-Dataset_1/m1/lobule4/0-Actin_Dendra_LD-1_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 4_Merged-2.tif',
               '../test-Dataset_1/m1/lobule4/1-Actin_Dendra_LD-2_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 4_Merged-2.tif',
               '../test-Dataset_1/m1/lobule4/2-Actin_Dendra_LD-3_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 4_Merged-2.tif',
               '../test-Dataset_1/m3/lobule2/0-Actin_Dendra_LD-1_merge_WD_M3_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
               '../test-Dataset_1/m3/lobule3/0-Actin_Dendra_LD-1_merge_WD_M3_12pm_dendra2_actin647_lipitox594.lif - TileScan 3_Merged-2.tif',
               '../test-Dataset_1/m3/lobule3/1-Actin_Dendra_LD-2_merge_WD_M3_12pm_dendra2_actin647_lipitox594.lif - TileScan 3_Merged-2.tif',
               '../test-Dataset_1/m3/lobule3/2-Actin_Dendra_LD-3_merge_WD_M3_12pm_dendra2_actin647_lipitox594.lif - TileScan 3_Merged-2.tif',
              
               '../test-Dataset_2/m1/lobule1/actin_dendra2_LD/0-Actin_Dendra_LD-1_M1_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
               '../test-Dataset_2/m1/lobule1/actin_dendra2_LD/1-Actin_Dendra_LD-2_M1_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
               '../test-Dataset_2/m1/lobule1/actin_dendra2_LD/2-Actin_Dendra_LD-3_M1_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-3.tif',
               '../test-Dataset_2/m1/lobule2/actin_dendra2_LD/0-Actin_Dendra_LD-1_M1_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 3_Merged-1.tif',
               '../test-Dataset_2/m1/lobule2/actin_dendra2_LD/1-Actin_Dendra_LD-2_M1_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 3_Merged-2.tif',
               '../test-Dataset_2/m1/lobule2/actin_dendra2_LD/2-Actin_Dendra_LD-3_M1_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 3_Merged-2.tif',
               '../test-Dataset_2/m2/lobule2/actin_dendra2_LD/0-Actin_Dendra_LD-1_M2_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
               '../test-Dataset_2/m2/lobule2/actin_dendra2_LD/1-Actin_Dendra_LD-2_M2_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
               '../test-Dataset_2/m2/lobule2/actin_dendra2_LD/2-Actin_Dendra_LD-3_M2_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
               '../test-Dataset_2/m3/lobule1/actin_dendra2_LD/0-Actin_Dendra_LD-1_M3_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
               '../test-Dataset_2/m3/lobule1/actin_dendra2_LD/1-Actin_Dendra_LD-2_M3_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
              '../test-Dataset_2/m3/lobule1/actin_dendra2_LD/2-Actin_Dendra_LD-3_M3_CNT_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
                
              '../test-Dataset_3/m1/lobule1/actin_dendra2_LD/2-Actin_Dendra_LD-3_M1_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-3.tif',
              '../test-Dataset_3/m1/lobule2/actin_dendra2_LD/0-Actin_Dendra_LD-1_M1_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
              '../test-Dataset_3/m2/lobule1/actin_dendra2_LD/0-Actin_Dendra_LD-1_M2_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-1.tif',
              '../test-Dataset_3/m2/lobule1/actin_dendra2_LD/1-Actin_Dendra_LD-2_M2_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
              '../test-Dataset_3/m2/lobule2/actin_dendra2_LD/0-Actin_Dendra_LD-1_M2_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
              '../test-Dataset_3/m3/lobule1/actin_dendra2_LD/0-Actin_Dendra_LD-1_M3_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
              '../test-Dataset_3/m3/lobule1/actin_dendra2_LD/1-Actin_Dendra_LD-2_M3_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
              '../test-Dataset_3/m3/lobule1/actin_dendra2_LD/2-Actin_Dendra_LD-3_M3_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 1_Merged-2.tif',
              '../test-Dataset_3/m3/lobule2/actin_dendra2_LD/0-Actin_Dendra_LD-1_M3_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
              '../test-Dataset_3/m3/lobule2/actin_dendra2_LD/1-Actin_Dendra_LD-2_M3_EF_EF_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif',
              
              ]
    
    lipid_model = models.CellposeModel(
        pretrained_model='models/lipid_droplet_model')
    mito_model = models.CellposeModel(
        pretrained_model='models/mito_model')
    cell_model = models.CellposeModel(
        pretrained_model='models/cell_model')
    
    start = time.time()
    for save_path, image_path in zip(save_paths, images):
        # segment.process_cells(image_path, cell_model, save_path)
        # segment.process_mito(image_path, mito_model, save_path)
        # segment.process_lipid_droplets(image_path, lipid_model, save_path)
        # segment.central_distance(save_path)
        # segment.portal_distance(save_path)
        # segment.cell_edge_distance(save_path)
        
        a = ex_prop.properties(save_path, 14.4042)
        a.cell_properties()
        visualize.label_ascini_cells(save_path)
        
        print(save_path)
        print(f'ran in {time.time() - start:.3f}')
        start = time.time()
        
        
        
        
        
        
        
        