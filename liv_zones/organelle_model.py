import numpy as np

from cellpose import models, utils
from cellpose.io import imread
import torch as torch

class OrganelleModel:
    def __init__(self, model_type):
        self.model_type = model_type

        if self.model_type == "cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="path_to_cell_model", device=torch.device("cuda")

            )
            self.flow = 0.4           #when control and fasted (STV) use 0.4
            self.cell_prob = 0.0
            self.diameter = 547.73     #when control and fasted (STV) use 254

        

        if self.model_type == "lipid_droplet_Small":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="path_to_ld_model", device=torch.device("cuda")
                
            )
            self.flow = 0.21
            self.cell_prob = 0.0
            self.diameter = 16  # also use 355 for large lipids
        
            
        
        if self.model_type == "lipid_droplet_Medium":
            self.pretrained_model = models.CellposeModel(pretrained_model="path_to_ld_model", device=torch.device("cuda")
                 
            )
            self.flow = 0.14
            self.cell_prob = 0.0
            self.diameter = 132  # also use 355 for large lipids
            
        if self.model_type == "lipid_droplet_Large":
            self.pretrained_model = models.CellposeModel(pretrained_model="path_to_ld_model", device=torch.device("cuda")
                 
            )
            self.flow = 0.4
            self.cell_prob = 0.0
            self.diameter = 645  # also use 355 for large lipids

        if self.model_type == "mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="path_to_mito_model", device=torch.device("cuda")
            )
            self.flow = 0.0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="path_to_peroxisome_model", device=torch.device("cuda")
            )
            self.flow = 0.0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "nuclei":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="path_to_nuclei_model", device=torch.device("cuda")
            )
            self.flow = 0.4
            self.cell_prob = 0.0
            self.diameter = None

    def segment(self, img_path, channel, save=True, save_path=""):

        image = imread(img_path)
        organelle_raw = image[channel, :, :]

        organelle_mask = self.pretrained_model.eval(
            organelle_raw,
            flow_threshold=self.flow,
            cellprob_threshold=self.cell_prob,
            diameter=self.diameter
        )

        if save is True:
            np.save(f"{save_path}/{self.model_type}_mask.npy", organelle_mask[0])

        return organelle_mask
