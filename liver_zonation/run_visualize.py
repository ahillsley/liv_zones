import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import visualize


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
    
    ds1_paths = save_paths[:9]
    ds2_paths = save_paths[9:21] 
    
    ds3_paths = save_paths[21:]
    
    
    a = visualize.group_ascini(save_paths, 'mitochondria_properties')
    a50 = np.percentile(a['aspect_ratio'], 50)
    a75 = np.percentile(a['aspect_ratio'], 75)
    a90 = np.percentile(a['aspect_ratio'], 90)
    
    plt.hist(a['aspect_ratio'], bins=500, range = (1,5))
    plt.vlines(a50, 0, 3000, color = 'k')
    plt.vlines(a75, 0, 3000, color = 'k')
    plt.vlines(a90, 0, 3000, color = 'k')
    plt.ylabel('Count')
    plt.xlabel('aspect ratio')
    plt.figure()
    
    
    data_set_1 = visualize.group_ascini(ds1_paths, 'average_properties_per_cell')
    data_set_2 = visualize.group_ascini(ds2_paths, 'average_properties_per_cell')
    data_set_3 = visualize.group_ascini(ds3_paths, 'average_properties_per_cell')
    
    aa
    trend_1 = visualize.get_trendline(data_set_1, window=200)
    trend_2 = visualize.get_trendline(data_set_2, window=200)
    trend_3 = visualize.get_trendline(data_set_3, window=200)
    
    visualize.plot_trends([trend_1, trend_2, trend_3])
    visualize.plot_individual_ascini_grouped(ds1_paths,
                                             ds2_paths,
                                             ds3_paths,
                                             color='g')
    visualize.plot_individual_ascini(ds1_paths, color='b')
    visualize.plot_individual_ascini(ds2_paths, color='orange')
    visualize.plot_individual_ascini(ds3_paths, color='g')
    
    
    visualize.single_cell_plots(
        image_path=images[24],
        save_path=save_paths[24],
        cell_num=50,
        mito_prop='area',
        ld_prop='area',
        cell_prop='mito_avg_area')
    
    
    visualize.plot_individual_cell(image_path=images[24],
                                   save_path=save_paths[24],
                                   cell_num=50)
    visualize.color_by_prop(image, prop_list, prop)
    
