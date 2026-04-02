import numpy as np
from pathlib import Path

from cellpose import models, utils
from cellpose.io import imread
import torch as torch

_MODELS_DIR = Path(__file__).parent / "models"


class OrganelleModel:
    def __init__(self, model_type):
        self.model_type = model_type

        if self.model_type == "cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=str(_MODELS_DIR / "cell_model"), device=torch.device("cuda")
            )
            self.flow = 0.4
            self.cell_prob = 0
            self.diameter = None

        if self.model_type == "lipid_droplet":
            # handels small lipid droplets
            self.pretrained_model = models.CellposeModel(
                pretrained_model=str(_MODELS_DIR / "lipid_droplet_model"),
                device=torch.device("cuda"),
            )
            self.flow = 0.3
            self.cell_prob = 0
            self.diameter = 23

        if self.model_type == "lipid_droplet_large":
            # handels large lipid droplets
            self.pretrained_model = models.CellposeModel(
                pretrained_model=str(_MODELS_DIR / "lipid_droplet_model"),
                device=torch.device("cuda"),
            )
            self.flow = 0.3
            self.cell_prob = 0
            self.diameter = 355

        if self.model_type == "mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=str(_MODELS_DIR / "mito_model"), device=torch.device("cuda")
            )
            self.flow = 0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=str(_MODELS_DIR / "peroxisome_model"), device=torch.device("cuda")
            )
            self.flow = 0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "nuclei":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=str(_MODELS_DIR / "nuclei_model"), device=torch.device("cuda")
            )
            self.flow = 0.4
            self.cell_prob = 0
            self.diameter = None

    def segment(self, img_path, channel=None, save=True, save_path=""):
        if isinstance(img_path, np.ndarray):
            organelle_raw = img_path
        else:
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
