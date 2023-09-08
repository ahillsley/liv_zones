import numpy as np

from cellpose import models, utils
from cellpose.io import imread


class Organelle_Model:
    def __init__(self, model_type):
        self.model_type = model_type

        if self.model_type == "cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/cell_model"
            )

        if self.model_type == "lipid_droplet":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/lipid_droplet_model"
            )

        if self.model_type == "mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/mito_model"
            )
        
        if self.model_type == "peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model='models/mito_model')

    def segment(self, img_path, channel, save=True, save_path=""):

        image = imread(img_path)
        organelle_raw = image[channel, :, :]

        organelle_mask = self.pretrained_model.eval(
            organelle_raw, flow_threshold=0.4, cellprob_threshold=0
        )
        if save is True:
            np.save(f"{save_path}/{self.model_type}_mask.npy", organelle_mask[0])

        return organelle_mask
