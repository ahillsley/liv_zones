# liv_zones — Project Summary

## Overview

`liv_zones` is a Python package for extracting and analyzing features of hepatic cells and organelles from fluorescent liver tissue microscopy images (confocal TIF stacks). It uses deep learning-based segmentation (Cellpose) and morphological analysis to quantify subcellular structures in the context of liver lobule anatomy (proximity to central vs. portal veins).

---

## Analysis Pipeline

1. **Preprocessing** — Segments cells and organelles from multi-channel images using pretrained Cellpose models; computes distance transforms relative to cell edges and veins.
2. **Organelle Feature Extraction** — Extracts morphological properties (area, aspect ratio, solidity, boundary distance, vein proximity, etc.) for mitochondria, lipid droplets, peroxisomes, and nuclei.
3. **Cell-Level Aggregation** — Summarizes organelle statistics (density, mean area, type distributions) at the per-cell level.
4. **Spatial Analysis** — Maps features along the portal-to-central vein axis of the liver lobule (acinus).
5. **Output** — Produces CSV tables for downstream statistical analysis and optional visualizations.

---

## Module Descriptions

| Module | Role |
|---|---|
| `preprocess.py` | Orchestrates mask generation; runs segmentation models and computes distance maps |
| `organelle_model.py` | Wrapper around Cellpose models for segmenting cells, mitochondria, lipid droplets, peroxisomes, and nuclei |
| `organelle.py` | Loads segmentation masks and computes per-organelle morphological features |
| `cell.py` | Aggregates organelle features to the cell level; produces per-cell summary statistics |
| `ascini.py` | Loads, filters, and visualizes cell properties across an acinus |
| `distance_to_veins.py` | Segments central and portal veins from cell masks; generates vein distance transforms |
| `crop.py` | Rotates and crops whole-lobule images into individual acinus samples |
| `track_z_pos.py` | Tracks cell centroids across z-slices using the Motile library |
| `correlationMatrix.py` | Computes and visualizes pairwise R² correlations between cell-level features |
| `paths.py` | Stores file system paths for all experimental groups (sex × diet condition) |
| `veincoords.py` | Stores pixel coordinates of central and portal vein locations per sample |
| `point.py` / `util.py` | Geometric helper classes and functions (2D points, rotation, bounding boxes) |
| `run.py` | Configurable top-level entry point for running the full pipeline |

---

## Key Features Extracted

- **Organelle-level**: area, perimeter, centroid, aspect ratio, solidity, distance to cell edge, distance to veins, ascinus position (−1 = central vein, +1 = portal vein)
- **Mitochondria** classified into 3 morphological types by aspect ratio
- **Lipid droplets** classified into 3 size classes by area
- **Cell-level**: organelle density, mean area, type fractions, and other aggregated statistics

---

## Experimental Design

Data is organized by animal sex (male/female) and diet condition (CNT / STV / WD), across multiple livers, lobules, and acini. Vein coordinates are stored per sample to enable spatial normalization.

---

## Dependencies

- `cellpose` — deep learning segmentation
- `scikit-image` — morphological analysis
- `numpy`, `pandas` — data handling
- `matplotlib` — visualization
- `motile` / `networkx` — cell tracking across z-slices (optional)
