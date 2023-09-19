from liv_zones import preprocess as pre
from liv_zones import organelle as org
from liv_zones import cell as c
from liv_zones.ascini import plot_properties, plot_cell, plot_ascinus_annotated
import time

""" 
Define scale and file paths
"""

scale = 22.187  # pixels per micron


image_paths = [
    "../../images/asinus_0/acinus0_z5_full.tif"
]

save_paths = [
    "../../images/asinus_0"
    ]


"""
Options on what to run
"""

# Do you want to run the preprocessing?
run_preprocessing = False

# comment out any features that you don't want re-calculated
feature_list = [
    "cell_mask",
    "mito_mask",
    "lipid_mask",
    "peroxisome_mask",
    "nuclei_mask",
    "central_dist",
    #"portal_dist",
    "boundary_distance",
]

channels = {"actin": 0, "nuclei": 1, "mito": 2, "lipid": 3, "peroxi": 4}

# Do you want to extract individual organelle features?
organelle_features = True

# comment out any organelles you dont want re-calculated
organelle_list = ["mitochondria", "lipid_droplets", "peroxisomes", "nuclei"]

# Do you want to calculate average features per cell?
cell_features = False

# Do you want an image of the ascinus with each cell labeled?
plot_labeled_ascinus = False

# Do you want to plot per cell properties over the length of the ascinus
plot_props = False

# Do you want to visualize a specific cell?
show_individual_cell = False
# the cell number corresponding to the labeled ascinus
cell_number = 10


"""
----------------------------------------
Run (don't change anything below)
----------------------------------------
"""

if __name__ == "__main__":
    start = time.time()
    for image_path, save_path in zip(image_paths, save_paths):
        print("processing image:")
        print(f"        {image_path}")
        print("saving data to:")
        print(f"        {save_path}")
        print("-" * 30)

        # pre-processing
        if run_preprocessing is True:
            pre.preprocessing(image_path, save_path, channels, feature_list)

        else:
            pre.file_check(save_path)

        # extract individual organelle features
        if organelle_features is True:

            print("extracting individual organelle features")
            org.organelle_features(
                path=save_path, scale=scale, organelle_list=organelle_list,
            )

        # extract average features per cell
        if cell_features is True:

            print("averaging features per cell")
            c.cell_features(save_path, scale=scale)

        # visualizing the results
        if plot_labeled_ascinus is True:
            plot_ascinus_annotated(save_path)

        if plot_props is True:
            plot_properties(save_path)

        if show_individual_cell is True:
            plot_cell(image_path, save_path, cell_num=cell_number)
        
        print(f'ran in {time.time()-start:.3f}s')
