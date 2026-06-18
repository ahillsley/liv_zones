.. blinx documentation master file

What is ``liv_zones`` ?
=======================

``liv_zones`` is a **single-cell organelle profiling** package, developed for the paper
`Spatial Organellomics Maps Cell State Diversity and Metabolic Adaptation in Tissues <https://www.biorxiv.org/content/10.1101/2024.12.06.627285v2>`_

``liv_zones`` combines deep learning–based segmentation
(`Cellpose <https://www.cellpose.org/>`_) with morphological analysis to detect
organelles — such as mitochondria, lipid droplets, and peroxisomes — and
quantify their shape, abundance, and spatial distribution on a per-cell basis.

By measuring organelle features across many cells, ``liv_zones`` makes it possible
to study how subcellular structure varies between cells and across spatial axes
of a tissue. It was originally developed to profile hepatocytes along the
portal-to-central vein axis of the liver lobule (the acinus), enabling direct
measurement of **liver zonation** at the organelle level, but the underlying
segmentation and feature-extraction workflow is general and can be extended to
other organelles, cell types, and tissues.

.. jupyter-execute::
  :hide-code:

  import liv_zones



Full Documentation:
===================

.. toctree::
  :maxdepth: 2

  install
  features
  tutorial
  api
