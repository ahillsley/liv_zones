import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cellpose.io import imread


def plot_individual_ascini(save_paths):
    temp = pd.read_csv(f'{save_paths[0]}average_properties_per_cell.csv')
    props = np.asarray(temp.keys().astype('str'))
    #props = np.delete(props, (0,1,2,3,4, 14))
    
    for prop in props:
        for path in save_paths:
            b = pd.read_csv(f'{path}average_properties_per_cell.csv')
            b = b[b['area'] > 200 ]
            b = b.sort_values('ascini_position')
            trend = b.rolling(50, min_periods=20).mean()
            plt.plot(trend['ascini_position'], trend[prop], c='g', linewidth=2)
        plt.title(prop, size=15)
        plt.xlabel('relative ascini position', size=15)
        plt.ylabel(prop, size=15)
        plt.figure()
            
    return


def group_ascini(paths, file_name):
    data_set = pd.read_csv(f'{paths[0]}{file_name}.csv')
    for path in paths[1:]:
        b = pd.read_csv(f'{path}{file_name}.csv')
        data_set = pd.concat((data_set, b))
        
    return data_set

def get_trendline(data_set, window):
    data_set = data_set[data_set['area'] > 200 ]
    data_set = data_set.sort_values('ascini_position')
   
    data_set_trend = data_set.rolling(window, min_periods=100).mean()
   
    return data_set_trend

def plot_trends(data_sets):
    colors = ['b', 'orange', 'g']
    props = np.asarray(data_sets[0].keys().astype('str'))
    props = np.delete(props, (0,1,2,3,4, 14))
    for prop in props:
        for i in range(len(data_sets)):
            plt.plot(data_sets[i]['ascini_position'], data_sets[i][prop],
                     c=colors[i], linewidth=5)
        plt.title(prop, size=15)
        plt.xlabel('ascini position', size=15)
        plt.ylabel(prop, size=15)
        plt.xticks(size=10)
        plt.yticks(size=10)
        plt.legend(['set_1', 'set_2', 'set_3'])
        plt.figure()
        
    return

def bbox(img, cell_num):
    a = np.where(img == cell_num)
    bbox = np.min(a[0]), np.max(a[0]), np.min(a[1]), np.max(a[1])
    return bbox

def plot_individual_cell(image_path, save_path, cell_num):
    image = imread(image_path)
    cell_mask = np.load(f'{save_path}cell_mask.npy')
    mito_mask = np.load(f'{save_path}mito_mask.npy')
    ld_mask = np.load(f'{save_path}ld_mask.npy')
    
    bound_box = bbox(cell_mask, cell_num)
    
    roi_cell = cell_mask[bound_box[0]:bound_box[1],bound_box[2]:bound_box[3]]
    roi_mito = mito_mask[bound_box[0]:bound_box[1],bound_box[2]:bound_box[3]]
    roi_ld = ld_mask[bound_box[0]:bound_box[1],bound_box[2]:bound_box[3]]
    roi_image = image[:,bound_box[0]:bound_box[1],bound_box[2]:bound_box[3]]
    
    fig, axs = plt.subplots(2, 3)
    axs[0,0].imshow(roi_cell)
    axs[0,1].imshow(roi_mito)
    axs[0,2].imshow(roi_ld)
    axs[1,0].imshow(roi_image[0,:,:])
    axs[1,1].imshow(roi_image[1,:,:])
    axs[1,2].imshow(roi_image[2,:,:])
    [axi.set_axis_off() for axi in axs.ravel()]
    plt.figure()
    
    return roi_cell, roi_mito, roi_ld, roi_image
    #return bound_box

def individual_cell_scatter(save_path, cell_num, prop):
    cell_props = pd.read_csv(f'{save_path}average_properties_per_cell.csv')
    
    b = cell_props.sort_values('ascini_position')
    trend = b.rolling(50, min_periods=20).mean()
    plt.plot(trend['ascini_position'], trend[prop], c='k', linewidth=2)
    plt.title(prop, size=15)
    plt.xlabel('relative ascini position', size=15)
    plt.ylabel(prop, size=15)
    
    plt.scatter(cell_props['ascini_position'], cell_props[prop])
    plt.scatter(cell_props['ascini_position'][cell_num - 1],
                cell_props[prop][cell_num - 1], c='orange')
    plt.ylabel(prop)
    plt.xlabel('Ascini position')
    plt.savefig(f'{save_path}single_cell_scatter.png', bbox_inches='tight')
    plt.show()
    return

def color_by_prop(image, prop_list, prop):
    new_img = np.zeros((image.shape[0], image.shape[1]))
    for i in np.unique(image)[1:]:
        prop_value = prop_list[prop][prop_list['label'] == i]
        new_img[image == i] = float(prop_value)
    
    return new_img

def single_cell_plots(image_path, save_path, cell_num, mito_prop, ld_prop, cell_prop):
    mito_props = pd.read_csv(f'{save_path}mitochondria_properties.csv')
    ld_props = pd.read_csv(f'{save_path}lipid_dropplet_properties.csv')
    
    roi_cell, roi_mito, roi_ld, roi_image = plot_individual_cell(
        image_path, save_path, cell_num)
    
    mito_colored = color_by_prop(roi_mito, mito_props, mito_prop)
    ld_colored = color_by_prop(roi_ld, ld_props, ld_prop)
    
    cmap = plt.cm.get_cmap("viridis").copy()
    cmap.set_under('#8D98A7')
    
    fig, axs = plt.subplots(2, 3)
    axs[0,0].imshow(roi_cell, cmap=cmap, vmin=1)
    axs[0,1].imshow(mito_colored, cmap=cmap, vmin=1)
    axs[0,2].imshow(ld_colored, cmap=cmap, vmin=1)
    axs[1,0].imshow(roi_image[0,:,:])
    axs[1,1].imshow(roi_image[1,:,:])
    axs[1,2].imshow(roi_image[2,:,:])
    [axi.set_axis_off() for axi in axs.ravel()]
    plt.savefig(f'{save_path}single_cell_subplots.png', bbox_inches='tight', dpi=400)
    
    plt.figure()
    individual_cell_scatter(save_path, cell_num, cell_prop)
    label_ascini_cells(save_path)
        
    return
    
    
def label_ascini_cells(save_path):
    cell_mask = np.load(f'{save_path}cell_mask.npy')
    cell_props = pd.read_csv(f'{save_path}average_properties_per_cell.csv')
    plt.imshow(cell_mask)
    plt.axis('off')
    for cell in cell_props['label']:
        plt.text(x=cell_props['centroid-1'][cell-1], y=cell_props['centroid-0'][cell-1], s=f'{cell}', size=2.5,
                 horizontalalignment='center', fontweight= 'semibold')
    plt.savefig(f'{save_path}labled_cells.png', dpi=300, bbox_inches='tight')
    
    return

# 
    


