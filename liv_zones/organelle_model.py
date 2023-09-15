import numpy as np

from cellpose import models, utils
from cellpose.io import imread


class OrganelleModel:
    def __init__(self, model_type):
        self.model_type = model_type

        if self.model_type == "cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/cell_model"
            )
            self.flow = 0.4
            self.cell_prob = 0
            self.diameter = None

        if self.model_type == "lipid_droplet":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/lipid_droplet_model"
            )
            self.flow = 0.3
            self.cell_prob = 0
            self.diameter = 23  # also use 355 for large lipids

        if self.model_type == "mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/mito_model"
            )
            self.flow = 0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/mito_model"
            )
            self.flow = 0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model.type == "nuclei":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/nuclei_model"
            )
            self.flow = 0.4
            self.cell_prob = 0
            self.diameter = None

    def segment(self, img_path, channel, save=True, save_path=""):

        image = imread(img_path)
        organelle_raw = image[channel, :, :]

        organelle_mask = self.pretrained_model.eval(
            organelle_raw,
            flow_threshold=self.flow,
            cellprob_threshold=self.cell_prob,
            diameter=self.diameter,
        )

        if save is True:
            np.save(f"{save_path}/{self.model_type}_mask.npy", organelle_mask[0])

        return organelle_mask
