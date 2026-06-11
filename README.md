# liv_zones

`liv_zones` is a **single-cell organelle profiling** package, developed for the paper 
[Spatial Organellomics Maps Cell State Diversity and Metabolic Adaptation in Tissues](https://www.biorxiv.org/content/10.1101/2024.12.06.627285v2)

`liv_zones` combines deep learning–based segmentation
([Cellpose](https://www.cellpose.org/)) with morphological analysis to detect
organelles — such as mitochondria, lipid droplets, and peroxisomes — and
quantify their shape, abundance, and spatial distribution on a per-cell basis.

By measuring organelle features across many cells, `liv_zones` makes it possible
to study how subcellular structure varies between cells and across spatial axes
of a tissue. It was originally developed to profile hepatocytes along the
portal-to-central vein axis of the liver lobule (the acinus), enabling direct
measurement of **liver zonation** at the organelle level, but the underlying
segmentation and feature-extraction workflow is general and can be extended to
other organelles, cell types, and tissues.

![Spatial Organellomics overview: tissue imaging by confocal microscopy is segmented at the tissue scale (single cells) and organelle scale (peroxisomes, mitochondria, lipid droplets) to yield single-cell organelle metrics — shape, size, amount, and position.](images/organelleomics_overview.png)

Full documentation is available
[here](https://ahillsley.github.io/liv_zones/).

## Installation

`liv_zones` requires Python and is installed from source. Using a fresh conda
environment is recommended:

```bash
conda create -n liv_zones python
conda activate liv_zones
git clone https://github.com/ahillsley/liv_zones.git
cd liv_zones
pip install .
```

## Sample Data

The raw imaging data from the study is available on figshare under a CC BY 4.0
license:

- [**Raw images dataset (Janelia figshare)**](https://janelia.figshare.com/articles/dataset/Raw_images_for_Spatial_single-cell_Organellomics_reveals_nutrient_dependent_hepatocyte_heterogeneity_and_predicts_pathophysiological_status_i_i_i_in_vivo_i_/31863250)
  — DOI [10.25378/janelia.31863250](https://doi.org/10.25378/janelia.31863250)

It contains raw confocal microscopy images of a representative **liver lobule**
and a **pancreas section** from a single control mouse, provided so you can test
the full segmentation → feature-extraction → spatial-analysis pipeline. Each
structure is labeled in its own fluorescence channel:

| Channel | Liver marker | Structure |
|:---:|:---|:---|
| C00 | mito-Dendra2 | Mitochondria |
| C01 | phalloidin | Actin / cell boundaries |
| C02 | LipidTOX | Lipid droplets |
| C03 | PMP70 | Peroxisomes |

The dataset is organized into three sets of multi-channel TIFs (one file per
z-slice per channel; ~234 files, **~88 GB** total):

| Folder | Contents | Z-slices | Channels | Per file | Total |
|:---|:---|:---:|:---:|:---:|:---:|
| `Con_liver1_lobule1_organelle_imgaing_dataset` | Full liver lobule z-stack | 0–50 (51) | 4 (C00–C03) | ~385 MB | ~79 GB |
| `small_subset_of_lobule1` | First 6 slices of that lobule, for quick testing | 0–5 (6) | 4 (C00–C03) | ~385 MB | ~9 GB |
| `male_pancreas2_region6` | Single pancreas plane (endocrine + exocrine) | 1 | 6 (ch00–ch05) | ~60 MB | ~0.4 GB |

The pancreas images carry the same 4 channels in addition to two extra channels labeling **insulin** and
**glucagon** to mark the endocrine region.


## Pipeline overview

The pipeline runs in three stages, mirrored by the three example notebooks
(`01_preprocessing`, `02_feature_extraction`, `03_visualization`).

### Input data

A single multi-channel image of a tissue section — a 2D plane (or z-slice)
stored as a TIF with one channel per labeled structure (e.g. actin/cell
boundaries, mitochondria, lipid droplets, peroxisomes). A `channels` dict maps
each structure to its channel index, and `scale` gives the image resolution in
pixels per micron. For spatial (zonation) analysis you also supply the pixel
coordinates of the two reference points that define the tissue axis — for liver,
the **central** and **portal** vein centers.

Each stage has **general** steps (organelle profiling on any tissue) and
**liver-specific** steps (spatial analysis along the portal-to-central vein
axis). The liver-specific steps can be swapped out or skipped when applying the
pipeline to other tissues.

### 1. Preprocessing & segmentation

Segments the cells and each organelle, then computes the distance maps that
encode where each pixel sits within the tissue. Orchestrated by
[`preprocess.preprocessing()`](liv_zones/preprocess.py).

*General — organelle segmentation:*
- [`organelle_model.OrganelleModel`](liv_zones/organelle_model.py) — wraps the per-organelle Cellpose models (one trained model per structure)
- [`preprocess.cell_edge_distance()`](liv_zones/preprocess.py) — distance from each pixel to the nearest cell boundary (within-cell position)

*Liver-specific — acinus geometry:*
- [`crop.dispCrop()`](liv_zones/crop.py) — crops a single acinus out of a whole-lobule image using the vein coordinates
- [`distance_to_veins.main()`](liv_zones/distance_to_veins.py) — distance maps from the central and portal veins

### 2. Feature extraction

Measures morphological features of every segmented organelle, then aggregates
them into per-cell summaries.

*General — organelle & cell features:*
- [`organelle.organelle_features()`](liv_zones/organelle.py) — per-organelle features (area, aspect ratio, solidity, …), one CSV per organelle type
- [`cell.cell_features()`](liv_zones/cell.py) — aggregates organelles to one row per cell (density, mean area, type fractions, …)

*Liver-specific — spatial position:*
- [`organelle.ascini_position()`](liv_zones/organelle.py) — each object's normalized position along the portal-to-central axis (−1 → +1)

### 3. Visualization & analysis

Explores how organelle properties vary across cells and along the tissue axis.

*General — per-cell & correlation views:*
- [`ascini.plot_cell()`](liv_zones/ascini.py) — detailed multi-channel view of an individual cell
- [`correlationMatrix.py`](liv_zones/correlationMatrix.py) — pairwise R² correlation heatmap between cell-level features

*Liver-specific — zonation trends:*
- [`ascini.plot_ascinus_annotated()`](liv_zones/ascini.py) — cell-segmentation map with each cell labeled by ID
- [`ascini.plot_properties()`](liv_zones/ascini.py) — per-cell feature trends along the acinus

### Outputs

All stages read from and write to a single output folder.

*General — organelle profiling:*

| Output | Produced by | Contents |
|:---|:---|:---|
| `*_mask.npy` | Stage 1 | Instance segmentation masks (cells, mitochondria, lipid droplets, peroxisomes) |
| `boundary_dist.npy` | Stage 1 | Distance from each pixel to the nearest cell boundary |
| `mitochondria_properties.csv`, `lipid_droplet_properties.csv`, `peroxisome_properties.csv` | Stage 2 | One row per organelle |
| `average_properties_per_cell.csv` | Stage 2 | One row per cell, for downstream statistics |
| correlation heatmap | Stage 3 | Pairwise R² between cell-level features |

*Liver-specific — acinus / zonation:*

| Output | Produced by | Contents |
|:---|:---|:---|
| `central_dist.npy`, `portal_dist.npy` | Stage 1 | Distance maps from the central and portal veins |
| `labeled_cells.png` | Stage 3 | Cell-segmentation map with each cell labeled by ID |
| `ascini_trends/` | Stage 3 | Per-cell feature trends along the acinus |


## Quickstart

The pipeline is driven by a handful of high-level functions. Each stage reads
from and writes to a single `save_path` folder, so the stages can be run
independently as long as the previous stage's outputs are present.

```python
from liv_zones import preprocess, organelle, cell
from liv_zones.ascini import plot_properties, plot_cell, plot_ascinus_annotated

image_path = "path/to/image.tif"   # multi-channel TIF
save_path  = "path/to/output"      # results are written here
scale      = 14.4024               # pixels per micron

# 1. Preprocessing & segmentation
#    Segments cells and organelles with Cellpose and computes the
#    vein / cell-boundary distance maps. Writes *_mask.npy and *_dist.npy.
preprocess.preprocessing(
    image_path,
    save_path,
    channels={"actin": 0, "mito": 1, "lipid": 2},
)

# 2. Feature extraction
#    organelle_features  -> one CSV per organelle (one row per organelle)
#    cell_features       -> average_properties_per_cell.csv (one row per cell)
organelle.organelle_features(
    save_path,
    scale,
    organelle_list=["mitos", "lipid_droplets", "peroxisomes"],
)
cell.cell_features(save_path, scale)

# 3. Visualization
plot_ascinus_annotated(save_path)        # cell map labeled by ID
plot_properties(save_path)               # feature trends along the acinus
plot_cell(save_path, cell_number=10)     # detailed view of one cell
```

| Function | Stage | What it does |
|:---|:---|:---|
| `preprocess.preprocessing()` | 1 | Run Cellpose segmentation and compute distance maps |
| `organelle.organelle_features()` | 2 | Per-organelle morphological features → CSV |
| `cell.cell_features()` | 2 | Aggregate organelle stats to one row per cell → CSV |
| `ascini.plot_*()` | 3 | Annotated cell maps, spatial trends, and per-cell views |

For batch processing, `run.py` in the package root wraps these calls with
configurable flags and loops over a list of images. The
[full tutorial](https://ahillsley.github.io/liv_zones/) walks through every
option.

## Tutorial

`notebooks` contains a series of interactive
Jupyter notebooks in written for biologists and new
Python users, with step-by-step explanations and inline visualizations. They
run on the small sample dataset shipped in `notebooks/sample_data/` — no data
download required.

```bash
conda activate liv_zones
cd notebooks/
jupyter notebook
```

Then work through the notebooks in order:

| Notebook | What it covers |
|:---|:---|
| `00_download_data.ipynb` | Locate and inspect a raw dataset (channels, z-slices) |
| `01_preprocessing.ipynb` | Load images, run segmentation, compute distance maps |
| `02_feature_extraction.ipynb` | Extract per-organelle and per-cell features |
| `03_visualization.ipynb` | Spatial trends, per-cell views, correlation heatmap |

Because the sample masks are pre-computed, you can start at
`02_feature_extraction.ipynb` even without a GPU.


## Citation

If you find this work useful for your research, please cite our
[paper](https://www.biorxiv.org/content/10.1101/2024.12.06.627285v2):

> Adhikari, R., Hillsley, A., Johnson, A.D., Gao, S.M., Espinosa-Medina, I.,
> Funke, J. and Feliciano, D., 2026. Spatial Organellomics Maps Cell State
> Diversity and Metabolic Adaptation in Tissues. _bioRxiv_ 2026

BibTeX:
```bibtex
@article{scOrganelleomics,
  author  = { Adhikari, Raghabendra and Hillsley, Alexander and Johnson, Alana D. and Gao, Shihong M. and Espinosa-Medina, Isabel and Funke, Jan and Feliciano, Daniel},
  title   = {Spatial Organellomics Maps Cell State Diversity and Metabolic Adaptation in Tissues},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.1101/2024.12.06.627285},
  url     = {https://www.biorxiv.org/content/10.1101/2024.12.06.627285v2}
}
```

## License

MIT
