import numpy as np

from cellpose import models, utils
from cellpose.io import imread
import torch as torch


class OrganelleModel:
    def __init__(self, model_type):
        self.model_type = model_type

        if self.model_type == "cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/cell_model", device=torch.device("cuda")
            )
            self.flow = 0.4
            self.cell_prob = 0
            self.diameter = None

        if self.model_type == "lipid_droplet":
            # handels small lipid droplets
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/lipid_droplet_model",
                device=torch.device("cuda"),
            )
            self.flow = 0.3
            self.cell_prob = 0
            self.diameter = 23

        if self.model_type == "lipid_droplet_large":
            # handels large lipid droplets
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/lipid_droplet_model",
                device=torch.device("cuda"),
            )
            self.flow = 0.3
            self.cell_prob = 0
            self.diameter = 355

        if self.model_type == "mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/mito_model", device=torch.device("cuda")
            )
            self.flow = 0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/peroxisome_model", device=torch.device("cuda")
            )
            self.flow = 0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "nuclei":
            self.pretrained_model = models.CellposeModel(
                pretrained_model="models/nuclei_model", device=torch.device("cuda")
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
