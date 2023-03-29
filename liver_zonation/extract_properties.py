from cellpose.io import imread
from skimage.measure import regionprops_table
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None



class properties:
    
    def __init__(self, save_path, scale=1):
        self.save_path = save_path
        self.cell_mask = np.load(f'{save_path}cell_mask.npy')
        self.cv_distance = np.load(f'{save_path}cv_distance.npy')
        self.pv_distance = np.load(f'{save_path}pv_distance.npy')
        self.cell_edge_distance = np.load(f'{save_path}cell_boundry_distance.npy')
        self.scale = scale
        
    def map_to_cell(self, props, mask):
        y_cords = np.asarray(props['centroid-0'], dtype='int')
        x_cords = np.asarray(props['centroid-1'], dtype='int')
        
        return mask[y_cords, x_cords]
    
    def mito_properties(self, min_mito_size=0):
        mito_mask = np.load(f'{self.save_path}mito_mask.npy')
        
        mito_props_dict = regionprops_table(mito_mask, properties =(
            'label',
            'area',
            'perimeter',
            'centroid',
            'axis_major_length',
            'axis_minor_length',
            'solidity'))
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
        
        mito_props['aspect_type_1'], \
        mito_props['aspect_type_2'], \
        mito_props['aspect_type_3'] = self.split_types(
            mito_props, 'aspect_ratio', (1.2, 2))
        
        mito_props.to_csv(f'{self.save_path}mitochondria_properties.csv')
        
        return mito_props
    
    def ascini_position(self, props):
        cv_distance = self.map_to_cell(props, self.cv_distance)
        pv_distance = self.map_to_cell(props, self.pv_distance)
        ascini_position = (pv_distance - cv_distance) / \
            (pv_distance + cv_distance)
        
        return ascini_position
    
    def lipid_droplet_properties(self, min_ld_size=0):
        ld_mask = np.load(f'{self.save_path}ld_mask.npy')
        
        ld_props_dict = regionprops_table(ld_mask, properties =(
            'label',
            'area',
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
        
        ld_props['aspect_type_1'], \
        ld_props['aspect_type_2'], \
        ld_props['aspect_type_3'] = self.split_types(
            ld_props, 'area', (500, 2000))
        
        ld_props.to_csv(f'{self.save_path}lipid_dropplet_properties.csv')
        
        return ld_props
    
    def split_types(self, organelle_list, prop, bounds):
        '''
        Split the given dataframe into types based on "prop" and cutoffs
        '''
        
        type_1 = organelle_list[prop] < bounds[0]
        type_3 = organelle_list[prop] >= bounds[1]
        type_2 = np.logical_or(type_1, type_3) == False
        
        return type_1, type_2, type_3
    
    def cell_properties(self):
        cell_props_dict = regionprops_table(
            self.cell_mask, properties = ('area', 'centroid', 'label'))
        cell_props_1 = pd.DataFrame.from_dict(cell_props_dict)
        cell_props_2 = pd.DataFrame(columns=['mito_density',
                                             'mito_avg_area',
                                             'mito_percent_total_area',
                                             'mito_distance_from_edge',
                                             'mito_aspect_ratio',
                                             'mito_solidity',
                                             'ld_density',
                                             'ld_avg_area',
                                             'ld_percent_total_area',
                                             'ld_distance_from_edge',
                                             'type_1_mito_density',
                                             'type_2_mito_density',
                                             'type_3_mito_density',
                                             'percent_type_1_mito',
                                             'percent_type_2_mito',
                                             'percent_type_3_mito',
                                             'type_1_ld_density',
                                             'type_2_ld_density',
                                             'type_3_ld_density',
                                             'percent_type_1_ld',
                                             'percent_type_2_ld',
                                             'percent_type_3_ld',
                                             'type_1_mito_dist_from_edge',
                                             'type_2_mito_dist_from_edge',
                                             'type_3_mito_dist_from_edge',
                                             'type_1_ld_dist_from_edge',
                                             'type_2_ld_dist_from_edge',
                                             'type_3_ld_dist_from_edge',
                                             ],
                                  index = range(0,len(cell_props_dict['label'])))
    
        cell_props = pd.concat((cell_props_1, cell_props_2), axis=1)
        cell_props['area'] = cell_props['area'] / self.scale**2
        
        mito_props = self.mito_properties()
        ld_props = self.lipid_droplet_properties()
        
        for cell in cell_props['label']:
            single_cell = mito_props[mito_props['cell_id']==cell]
            
            # Standard properties
            # ----------------------------------------------------
            cell_props['mito_density'][cell-1] = len(single_cell) / \
                (cell_props['area'][cell-1] / self.scale**2)
                
            cell_props['mito_avg_area'][cell-1] = np.mean(
                single_cell['area'] / self.scale**2)
        
            cell_props['mito_percent_total_area'][cell-1] = np.sum(
                single_cell['area']) / cell_props['area'][cell-1]
            
            cell_props['mito_aspect_ratio'][cell-1] = np.mean(
                single_cell['aspect_ratio'])
            
            cell_props['mito_solidity'][cell-1] = np.mean(
                single_cell['solidity'])
        
            single_cell_ld = ld_props[ld_props['cell_id'] ==cell]
            
            cell_props['ld_density'][cell-1] = len(single_cell_ld) / \
                (cell_props['area'][cell-1] / self.scale**2)
                
            cell_props['ld_avg_area'][cell-1] = np.mean(
                single_cell_ld['area'] / self.scale**2)
            
            cell_props['ld_percent_total_area'][cell-1] = np.sum(
                single_cell_ld['area']) / cell_props['area'][cell-1]
            
            # Position within cell and ascini
            # ----------------------------------------------------
            cell_props['mito_distance_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1], single_cell)
            
            cell_props['ld_distance_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1], single_cell_ld)
            
            cell_props['ascini_position'] = self.ascini_position(cell_props)
            
            # Specific organelle type properties
            # ----------------------------------------------------
            cell_props['type_1_mito_density'][cell-1] = np.sum(
                single_cell['aspect_type_1']) / (
                    cell_props['area'][cell-1] / self.scale**2)
            
            cell_props['type_2_mito_density'][cell-1] = np.sum(
                single_cell['aspect_type_2']) / (
                    cell_props['area'][cell-1] / self.scale**2)
            
            cell_props['type_3_mito_density'][cell-1] = np.sum(
                single_cell['aspect_type_3']) / (
                    cell_props['area'][cell-1] / self.scale**2)
            
            cell_props['percent_type_1_mito'][cell-1] = np.sum(
                single_cell['aspect_type_1']) / len(single_cell)
            
            cell_props['percent_type_2_mito'][cell-1] = np.sum(
                single_cell['aspect_type_2']) / len(single_cell)
            
            cell_props['percent_type_3_mito'][cell-1] = np.sum(
                single_cell['aspect_type_3']) / len(single_cell)
            
            cell_props['type_1_mito_dist_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1],
                                    single_cell[single_cell['aspect_type_1']])
                
            cell_props['type_2_mito_dist_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1],
                                    single_cell[single_cell['aspect_type_2']])
            
            cell_props['type_3_mito_dist_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1],
                                    single_cell[single_cell['aspect_type_3']])
            
            # lipid_droplets
            cell_props['type_1_ld_density'][cell-1] = np.sum(
                single_cell_ld['aspect_type_1']) / (
                    cell_props['area'][cell-1] / self.scale**2)
            
            cell_props['type_2_ld_density'][cell-1] = np.sum(
                single_cell_ld['aspect_type_2']) / (
                    cell_props['area'][cell-1] / self.scale**2)
            
            cell_props['type_3_ld_density'][cell-1] = np.sum(
                single_cell_ld['aspect_type_3']) / (
                    cell_props['area'][cell-1] / self.scale**2)
            
            cell_props['percent_type_1_ld'][cell-1] = np.sum(
                single_cell_ld['aspect_type_1']) / len(single_cell_ld)
            
            cell_props['percent_type_2_ld'][cell-1] = np.sum(
                single_cell_ld['aspect_type_2']) / len(single_cell_ld)
            
            cell_props['percent_type_3_ld'][cell-1] = np.sum(
                single_cell_ld['aspect_type_3']) / len(single_cell_ld)
            
            cell_props['type_1_ld_dist_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1],
                                    single_cell_ld[
                                        single_cell_ld['aspect_type_1']])
                
            cell_props['type_2_ld_dist_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1],
                                    single_cell_ld[
                                        single_cell_ld['aspect_type_2']])
            
            cell_props['type_3_ld_dist_from_edge'][cell-1] = \
                self.dist_from_edge(cell_props.iloc[cell-1],
                                    single_cell_ld[
                                        single_cell_ld['aspect_type_3']])
                
        cell_props.to_csv(f'{self.save_path}average_properties_per_cell.csv')
        
        return cell_props
    
    def dist_from_edge(self, cell, organelle_list):
        cc_x = float(cell['centroid-0'])
        cc_y = float(cell['centroid-1'])
        
        o_s = np.asarray([organelle_list['centroid-0'],
                          organelle_list['centroid-1']])
        
        cc_s = np.asarray([np.repeat(cc_x, o_s.shape[1]),
                            np.repeat(cc_y, o_s.shape[1])])
        # in pixels
        dist_to_center = np.linalg.norm(cc_s - o_s, axis=0)
        dist_to_edge = np.asarray(organelle_list['boundry_dist'])
        rel_position = dist_to_center / (dist_to_center + dist_to_edge)
        
        return np.mean(rel_position)
        