from cellpose.io import imread
from skimage.measure import regionprops_table
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None



class properties:
    
    def __init__(self, save_path, scale=1):
        self.cell_mask = np.load(f'{save_path}_cell_mask.npy')
        self.cv_distance = np.load(f'{save_path}cv_distance.npy')
        self.pv_distance = np.load(f'{save_path}pv_distance.npy')
        self.cell_edge_distance = np.load(f'{save_path}_cell_boundry_distance.npy')
        self.scale = scale
        
    def map_to_cell(self, props, mask):
        y_cords = np.asarray(props['centroid-0'], dtype='int')
        x_cords = np.asarray(props['centroid-1'], dtype='int')
        
        return mask[y_cords, x_cords]
    
    def mito_properties(self, save_path, min_mito_size=0):
        mito_mask = np.load(f'{save_path}_mito_mask.npy')
        
        mito_props_dict = regionprops_table(mito_mask, properties =('area',
                                                        'perimeter',
                                                        'centroid',
                                                        'axis_major_length',
                                                        'axis_minor_length'))
        mito_props = pd.DataFrame.from_dict(mito_props_dict)
        mito_props['cell_id'] = self.map_to_cell(mito_props, self.cell_mask)
        mito_props['aspect_ratio'] = mito_props['axis_major_length'] / \
            mito_props['axis_minor_length']
        
        # use formula to calc boundry to cell edge
        mito_props['boundry_dist'] = self.map_to_cell(
            mito_props, self.cell_edge_distance)   
            
        
        # use new ascini distance measure
        # portal vein = 1, central vein = -1
        mito_props['ascini_position'] = self.ascini_position(mito_props)
        
        mito_props.to_csv(f'{save_path}mitochondria_properties.csv')
        
        return mito_props
    
    def ascini_position(self, props):
        cv_distance = self.map_to_cell(props, self.cv_distance)
        pv_distance = self.map_to_cell(props, self.pv_distance)
        ascini_position = (pv_distance - cv_distance) / \
            (pv_distance + cv_distance)
        
        return ascini_position
    
    def lipid_droplet_properties(self, save_path, min_ld_size=0):
        ld_mask = np.load(f'{save_path}ld_mask.npy')
        
        ld_props_dict = regionprops_table(ld_mask, properties =('area',
                                                        'perimeter',
                                                        'centroid',
                                                        'axis_major_length',
                                                        'axis_minor_length'))
        ld_props = pd.DataFrame.from_dict(ld_props_dict)
        ld_props['cell_id'] = self.map_to_cell(ld_props, self.cell_mask)
        ld_props['aspect_ratio'] = ld_props['axis_major_length'] / \
            ld_props['axis_minor_length']
            
        ld_props['boundry_dist'] = self.map_to_cell(
            ld_props, self.cell_edge_distance)  
        
        ld_props['ascini_position'] = self.ascini_position(ld_props)
        
        ld_props.to_csv(f'{save_path}lipid_dropplet_properties.csv')
        
        return ld_props
    
    def cell_properties(self, save_path):
        cell_props_dict = regionprops_table(
            self.cell_mask, properties = ('area', 'centroid', 'label'))
        cell_props_1 = pd.DataFrame.from_dict(cell_props_dict)
        cell_props_2 = pd.DataFrame(columns=['mito_density',
                                             'mito_avg_area',
                                             'mito_percent_total_area',
                                             'mito_distance_from_edge',
                                             'mito_aspect_ratio',
                                             'ld_density',
                                             'ld_avg_area',
                                             'ld_percent_total_area',
                                             'ld_distance_from_edge',
                                             ],
                                  index = range(0,len(cell_props_dict['label'])))
    
        cell_props = pd.concat((cell_props_1, cell_props_2), axis=1)
        cell_props['area'] = cell_props['area'] / self.scale**2
        
        mito_props = self.mito_properties(save_path)
        ld_props = self.lipid_droplet_properties(save_path)
        
        for cell in cell_props['label']:
            single_cell = mito_props[mito_props['cell_id']==cell]
            cell_props['mito_density'][cell-1] = len(single_cell) / \
                (cell_props['area'][cell-1] / self.scale**2)
                
            cell_props['mito_avg_area'][cell-1] = np.mean(
                single_cell['area'] / self.scale**2)
        
            cell_props['mito_percent_total_area'][cell-1] = np.sum(
                single_cell['area']) / cell_props['area'][cell-1]
            
            cell_props['mito_aspect_ratio'][cell-1] = np.mean(
                single_cell['aspect_ratio'])
        
            max_dist = np.max(self.cell_edge_distance[self.cell_mask==cell])
            cell_props['mito_distance_from_edge'][cell-1] = np.mean(
                (max_dist - single_cell['boundry_dist']) / \
                max_dist) / self.scale
            
            single_cell_ld = ld_props[ld_props['cell_id'] ==cell]
            
            cell_props['ld_density'][cell-1] = len(single_cell_ld) / \
                (cell_props['area'][cell-1] / self.scale**2)
                
            cell_props['ld_avg_area'][cell-1] = np.mean(
                single_cell_ld['area'] / self.scale**2)
            
            cell_props['ld_percent_total_area'][cell-1] = np.sum(
                single_cell_ld['area']) / cell_props['area'][cell-1]
            
            cell_props['ld_distance_from_edge'][cell-1] = np.mean(
                (max_dist - single_cell_ld['boundry_dist']) / \
                max_dist) / self.scale
            
            cell_props['ascini_position'] = self.ascini_position(cell_props)
        
        cell_props.to_csv(f'{save_path}average_properties_per_cell.csv')
        
        return cell_props