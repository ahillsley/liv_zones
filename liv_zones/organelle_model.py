import os

import numpy as np

from cellpose import models, utils
from cellpose.io import imread
import torch as torch
from skimage.segmentation import relabel_sequential

# Resolve the bundled models relative to this file so paths work regardless
# of the current working directory (e.g. when run from notebooks/).
_MODELS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")


def _model_path(name):
    return os.path.join(_MODELS_DIR, name)


class OrganelleModel:
    def __init__(self, model_type):
        self.model_type = model_type

        if self.model_type == "cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("cell_model"), device=torch.device("cuda")

            )
            self.flow = 0.4           #when control and fasted (STV) use 0.4
            self.cell_prob = 0.0
            self.diameter = 547.73     #when control and fasted (STV) use 254

        

        if self.model_type == "lipid_droplet_Small":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("lipid_droplet_model"), device=torch.device("cuda")

            )
            self.flow = 0.21
            self.cell_prob = 0.0
            self.diameter = 16  # also use 355 for large lipids
        
            
        
        if self.model_type == "lipid_droplet_Medium":
            self.pretrained_model = models.CellposeModel(pretrained_model=_model_path("lipid_droplet_model"), device=torch.device("cuda")

            )
            self.flow = 0.14
            self.cell_prob = 0.0
            self.diameter = 132  # also use 355 for large lipids
            
        if self.model_type == "lipid_droplet_Large":
            self.pretrained_model = models.CellposeModel(pretrained_model=_model_path("lipid_droplet_model"), device=torch.device("cuda")

            )
            self.flow = 0.4
            self.cell_prob = 0.0
            self.diameter = 645  # also use 355 for large lipids

        if self.model_type == "mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("mito_model"), device=torch.device("cuda")
            )
            self.flow = 0.0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("peroxisome_model"), device=torch.device("cuda")
            )
            self.flow = 0.0
            self.cell_prob = -0.1
            self.diameter = None

        if self.model_type == "nuclei":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("nuclei_model"), device=torch.device("cuda")
            )
            self.flow = 0.4
            self.cell_prob = 0.0
            self.diameter = None

        # ------------------------------------------------------------------
        # Pancreas models (see preprocess.preprocessing(tissue="pancreas"))
        # ------------------------------------------------------------------
        if self.model_type == "pancreas_cell":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("pancreas_cell_model"),
                device=torch.device("cuda"),
            )
            self.flow = 0.8
            self.cell_prob = 0.0
            self.diameter = 328

        if self.model_type == "pancreas_mito":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("pancreas_mito_model"),
                device=torch.device("cuda"),
            )
            self.flow = 0.0
            self.cell_prob = -0.15
            self.diameter = None

        if self.model_type == "pancreas_peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("pancreas_peroxisome_model"),
                device=torch.device("cuda"),
            )
            self.flow = 0.5
            self.cell_prob = 0.0
            self.diameter = 13.62

        if self.model_type == "pancreas_large_peroxisome":
            self.pretrained_model = models.CellposeModel(
                pretrained_model=_model_path("pancreas_large_peroxisome_model"),
                device=torch.device("cuda"),
            )
            self.flow = 0.5
            self.cell_prob = 0.0
            self.diameter = 42.76

    # ------------------------------------------------------------------
    # Tiled segmentation with IoU-based stitching, used for the large
    # pancreas cell images. Tiles the image, runs Cellpose per tile, and
    # stitches objects with greedy one-to-one IoU assignment.
    # ------------------------------------------------------------------
    def _segment_tiled_pancreas_cell(
        self,
        img2d,
        tile_size=2048,
        overlap=768,
        IoU_thresh=0.30,
        min_overlap_pixels=150,
        tile_frac_thresh=0.25,
    ):
        """
        Tile the large image, run Cellpose on each tile (crop-style), and stitch using
        IoU-based greedy one-to-one assignment with safe relabeling at the end.

        Parameters
        ----------
        img2d : 2D numpy uint8
        tile_size : int
            tile edge length
        overlap : int
            overlap in pixels between adjacent tiles
        IoU_thresh : float
            minimum IoU to consider a tile object matched to a global object
        min_overlap_pixels : int
            minimum number of overlapping pixels to consider candidate
        tile_frac_thresh : float
            minimum fraction of the tile-object area that must overlap the global object

        Returns
        -------
        full_mask : 2D int32 label image (sequential labels 1..N)
        """
        if img2d.ndim != 2:
            raise ValueError(f"Expected 2D image, got shape {img2d.shape}")

        H, W = img2d.shape
        step = max(1, tile_size - overlap)
        full_mask = np.zeros((H, W), dtype=np.int32)
        next_id = 1

        for y_start in range(0, H, step):
            y_end = min(y_start + tile_size, H)
            y_start = max(0, y_end - tile_size)

            for x_start in range(0, W, step):
                x_end = min(x_start + tile_size, W)
                x_start = max(0, x_end - tile_size)

                tile = img2d[y_start:y_end, x_start:x_end]

                # segment tile (crop-style). Use tile=False because we tile manually.
                try:
                    out = self.pretrained_model.eval(
                        tile,
                        flow_threshold=self.flow,
                        cellprob_threshold=self.cell_prob,
                        diameter=self.diameter,
                        tile=False,
                        normalize={"tile_norm": True},
                    )
                except TypeError:
                    out = self.pretrained_model.eval(
                        tile,
                        flow_threshold=self.flow,
                        cellprob_threshold=self.cell_prob,
                        diameter=self.diameter,
                        tile=False,
                        normalize=True,
                    )

                tile_mask = out[0].astype(np.int32)
                if tile_mask.max() == 0:
                    continue

                region = full_mask[y_start:y_end, x_start:x_end]

                # Tile labels (nonzero)
                tile_labels = np.unique(tile_mask)
                tile_labels = tile_labels[tile_labels != 0]
                if tile_labels.size == 0:
                    continue

                # Global labels present in the region
                global_labels = np.unique(region)
                global_labels = global_labels[global_labels != 0]

                # Precompute areas
                tile_areas = {lab: int((tile_mask == lab).sum()) for lab in tile_labels}
                global_areas = {gl: int((region == gl).sum()) for gl in global_labels}

                # Build candidate matches based on overlap, IoU, AND tile-fraction
                candidates = []  # (iou, tile_lab, global_lab, inter)
                for lab in tile_labels:
                    mask_lab = (tile_mask == lab)
                    t_area = tile_areas[lab]
                    if t_area == 0:
                        continue
                    overlapped = region[mask_lab]
                    overlapped = overlapped[overlapped != 0]
                    if overlapped.size == 0:
                        continue
                    glabs, counts = np.unique(overlapped, return_counts=True)
                    for g, cnt in zip(glabs, counts):
                        inter = int(cnt)
                        union = t_area + global_areas.get(int(g), 0) - inter
                        if union <= 0:
                            continue
                        iou = inter / union
                        frac = inter / t_area
                        # require absolute overlap AND IoU AND tile-fraction
                        if inter >= min_overlap_pixels and iou >= IoU_thresh and frac >= tile_frac_thresh:
                            candidates.append((iou, int(lab), int(g), inter))

                # Greedy one-to-one assignment by descending IoU
                candidates.sort(reverse=True, key=lambda x: x[0])
                mapping = {}  # tile_lab -> global_lab
                used_globals = set()
                for iou, t_lab, g_lab, inter in candidates:
                    if t_lab in mapping:
                        continue
                    if g_lab in used_globals:
                        continue
                    mapping[t_lab] = g_lab
                    used_globals.add(g_lab)

                # compute how many tile labels will become new global ids for this tile
                new_assigned = sum(1 for lab in tile_labels if lab not in mapping)

                # Debug print per tile
                print(f"tile {(y_start, x_start)}: candidates {len(candidates)}, assigned {len(mapping)}, new_ids {new_assigned}")

                # Create relabeled tile: mapped => existing global id; unmapped => new id
                relabeled = np.zeros_like(tile_mask, dtype=np.int32)
                for lab in tile_labels:
                    if lab in mapping:
                        relabeled[tile_mask == lab] = mapping[lab]
                    else:
                        relabeled[tile_mask == lab] = next_id
                        next_id += 1

                # Paste into full mask WITHOUT overwriting other non-zero labels:
                yy0, yy1 = y_start, y_end
                xx0, xx1 = x_start, x_end
                region = full_mask[yy0:yy1, xx0:xx1]

                # 1) Fill strictly empty pixels only
                write_mask = (region == 0) & (relabeled > 0)
                region[write_mask] = relabeled[write_mask]

                # 2) For mapped globals, allow relabeled pixels mapped to that global to extend label
                #    but ONLY where region is empty or already that same global label (do not overwrite other labels)
                for g_lab in used_globals:
                    extend_mask = (relabeled == g_lab) & ((region == 0) | (region == g_lab))
                    if np.any(extend_mask):
                        region[extend_mask] = g_lab

                full_mask[yy0:yy1, xx0:xx1] = region

        # Debug print before final sequential relabel
        print("global unique labels before relabel:", np.unique(full_mask).size)

        # Final cleanup: relabel existing labels sequentially without merging touching objects
        full_mask, _, _ = relabel_sequential(full_mask)
        full_mask = full_mask.astype(np.int32)

        return full_mask

    def segment(self, img_path, channel, save=True, save_path=""):

        image = imread(img_path)
        organelle_raw = image[channel, :, :]

        # compatibility: ensure diam_labels exists for different cellpose versions
        if not hasattr(self.pretrained_model, "diam_labels"):
            self.pretrained_model.diam_labels = None

        # Use tiled stitching only for the large pancreas cell images.
        if self.model_type == "pancreas_cell":
            mask2d = self._segment_tiled_pancreas_cell(
                organelle_raw,
                tile_size=2048,
                overlap=512,
                IoU_thresh=0.30,
                min_overlap_pixels=150,
                tile_frac_thresh=0.10,
            )
            organelle_mask = (mask2d, None, None)
        else:
            organelle_mask = self.pretrained_model.eval(
                organelle_raw,
                flow_threshold=self.flow,
                cellprob_threshold=self.cell_prob,
                diameter=self.diameter,
            )

        if save is True:
            # Pancreas masks are saved under the same canonical names as the
            # liver masks so all downstream modules stay tissue-agnostic.
            save_name = {
                "pancreas_cell": "cell",
                "pancreas_mito": "mito",
                "pancreas_peroxisome": "peroxisome",
                "pancreas_large_peroxisome": "large_peroxisome",
            }.get(self.model_type, self.model_type)
            np.save(f"{save_path}/{save_name}_mask.npy", organelle_mask[0])

        return organelle_mask
