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
   
  1. ``run_preprocessing``: will check to make sure that all the files needed for the analysis are present in the specific `save_path` folder, and will create any that are missing. The necessary files are:

    * ``cell_mask.npy``: an instance segmentation mask of all cells in the image

    * ``mito_mask.npy``: an instance segmentation mask of all mitochondria in the image

    * ``lipid_droplet_mask.npy``: an instance segmentation mask of all lipid droplets in the image

    * ``cv_distance.npy``: a distance transform from the central vein.

    * ``pv_distance.npy``: a distance transform from the portal vein

    * ``boundry_distance.npy``: a distance transform from the boundry of each cell

    * It is important that each of these files are named exactly as stated above, if not the program will not recognize them

  2. ``organelle_features``: will extract features 

