Running liv_zones
=================

Before running the liv_zones processing it is important to make sure that your images are in the correct format 
and that your folders are organized correctly

* Images are expected to be saved as .tif files, although we are working on also supporting .zarr 

* Images channels are expected to be organized as follows:

  0. Actin

  1. Mitochondria

  2. Lipid Droplets

  * additional channels are allowed, but they will not be processed

* There must be a seperate folder to store the results of each input image

  * outputs are saved with pre-defined names and will override previous outputs if that are any already in 
    the designated folder
  * I personally recommend creating a seperate folder for each image, then storing the original image as well
    as any outputs within

* The liv_zones analysis also requires 2 other inputs, saved in the output folder:

  * ``cv_distnace.tif``: a binary image segmentation of the central vein, where the central vein is marked by 0s

  * ``pv_distance.tif``: a binary image segmentation of the protal vein, where the portal vein marked by 0s

  * these are used to caclulate the realtive position of each organelle and cell along the ascinus

  * lastly, it is important that these are saved with these exact names, otherwise the program will not recognize them


In order to run the liv_zones processing, we use the ``run.py`` script

1. Open run.py in any text editor

.. jupyter-execute::
  
   from liv_zones import preprocess as pre
   from liv_zones import organelle as org
   from liv_zones import cell as c
   from liv_zones.ascini import plot_properties, plot_cell, plot_ascinus_annotated


2. Next we specify the scale and the needed file paths

   * ``scale``: the number of pixels per micron in the image, easiest to find using FIJI

   * ``image_paths``: a list of file paths to specific images for processing

   * ``save_paths``: a list of fie paths to the folders where results for each image are saved

    * the image_paths and save_paths lists are paired, ie. the results of the first image will be saved to the first
      save_path and so on... 

    * double check to make sure that both lists are the same length


.. jupyter-execute::

    scale =  14.4024 # pixels per micron 
    image_file_paths = [ 'path/to/first/image.tif',
                         'path/to/second/image.tif'
                         ]

    save_file_paths = [ 'path/to/first/folder',
                        'path/to/second/folder',
                        ]

3. Next we specify which analyses would like to run
   
  1. ``run_preprocessing``: `True` or `False`, If `True`, will check to make sure that all the files needed for the analysis are present in the specific `save_path` folder, and will create any that are missing. The necessary files are:

    * ``cell_mask.npy``: an instance segmentation mask of all cells in the image

    * ``mito_mask.npy``: an instance segmentation mask of all mitochondria in the image

    * ``lipid_droplet_mask.npy``: an instance segmentation mask of all lipid droplets in the image

    * ``cv_distance.npy``: a distance transform from the central vein.

    * ``pv_distance.npy``: a distance transform from the portal vein

    * ``boundry_distance.npy``: a distance transform from the boundry of each cell

    * It is important that each of these files are named exactly as stated above, if not the program will not recognize them

  .. jupyter-execute::

    # Do you want to run the preprocessing?
    run_preprocessing = False

    # comment out any features that you don't want re-calculated
    feature_list = [
      'cell_mask',
      'mito_mask',
      'lipid_droplet_mask',
      'cv_distance',
      'pv_distance',
      'boundry_distance'
      ]


  2. ``organelle_features``: `True` or `False`, If `True`, will extract features of the organelles provided in `organelle_list` and save them in a csv file

    * ``organelle_list``: possible options are `mitochondria` or `lipid_droplets`, comment out any that you do not wish to calculate

   .. jupyter-execute::

     #  Do you want to extract individual organelle features?
     organelle_features = False

     # comment out any organelles you dont want re-calculated
     organelle_list = [
        'mitochondria',
        'lipid_droplets',
        ]

  3. Visualization options:

    * ``plot_labled_ascinus``: `True` or `False`: If `True` save an image of the ascinus with each cell labeled with its ID #

    * ``plot_props``: `True` or `False`: If `True` generate plots of average properties per cell as a function of ascinus position. These plots will be saved in a sub-folder `ascini_trends`

    * ``show_individual_cell``: `True` or `False`: 

   .. jupyter-execute::

     # Do you want an image of the ascinus with each cell labeled?
     plot_labeled_ascinus = False

     # Do you want to plot per cell properties over the length of the ascinus
     plot_props = False

     # Do you want to visualize a specific cell?
     show_individual_cell = False

     # the cell number corresponding to the labeled ascinus
     cell_number = 10 

4. Run the analysis by calling saving your changes and running ``run.py`` from the terminal

  .. code-block:: bash

    python run.py
