import pandas as pd

import preprocess as pre
import organelle as org
import cell as c


""" 
Define constants and file paths
"""

scale = 14.4024  # pixels per micron

ld_area_split = (2.41, 9.64)

mito_aspect_split = (1.2, 2)


image_paths = [
    "//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv2/Lobule2/z05/acinus0_z5_full.tif"
]

save_paths = ["//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv2/Lobule2/z05"]


"""
Options
"""

run_preprocessing = True
organelle_features = False
cell_features = False


"""
----------------------------------------
Run (don't change anything below)
----------------------------------------
"""

if __name__ == "__main__":

    for image_path, save_path in zip(image_paths, save_paths):
        print("processing image:")
        print(f"        {image_path}")
        print("saving data to:")
        print(f"        {save_path}")
        print("-" * 30)

        # pre-processing
        if run_preprocessing is True:
            pre.preprocess(image_path, save_path)

        else:
            pre.file_check(save_path)

        # extract individual organelle features
        if organelle_features is True:

            print("extracting individual organelle features")
            org.organelle_features(
                path=save_path,
                scale=scale,
                mito_aspect_split=mito_aspect_split,
                ld_area_split=ld_area_split,
            )

        # extract average features per cell
        if cell_features is True:

            print("averaging features per cell")
            c.cell_features(save_path, scale=scale)
