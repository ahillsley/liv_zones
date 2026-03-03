# Script to write vein coordinate pickle files for each lobule.
#
# For each lobule directory, set:
#   CV_coords: [x, y] pixel coordinates of the central vein
#   PV_coords: [[x, y], ...] list of portal vein coordinates
#
# Then run this script once to write vein_coords.pickle into each lobule directory.
# These pickle files are read by the pipeline to compute zone distances.
#
# Repeat the os.chdir / CV_coords / PV_coords / pickle.dump block for each lobule.

import os
import pickle


# Example: one lobule
os.chdir('/path/to/your/data/Sex/Diet/Liv1/Lobule1')
CV_coords = [1000, 1000]   # replace with actual central vein pixel coordinates
PV_coords = [[500, 500], [1500, 500], [1000, 1500]]  # replace with actual portal vein coordinates
with open('vein_coords.pickle', 'wb') as f:
    pickle.dump([CV_coords, PV_coords], f, pickle.HIGHEST_PROTOCOL)

# Add more blocks below for each additional lobule, e.g.:
# os.chdir('/path/to/your/data/Sex/Diet/Liv1/Lobule2')
# CV_coords = [...]
# PV_coords = [...]
# with open('vein_coords.pickle', 'wb') as f:
#     pickle.dump([CV_coords, PV_coords], f, pickle.HIGHEST_PROTOCOL)
