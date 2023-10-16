import numpy as np
import glob
import warnings
from scipy.ndimage import distance_transform_edt  # uninstall and re-install

warnings.filterwarnings("ignore")

from cellpose import models, utils
from cellpose.io import imread


from liv_zones import organelle_model as org_models
from liv_zones.distance_to_veins import main as vein_dist


# Always define save path not including the last /

def preprocessing(image_path, save_path, channels, feature_list=None):
    # function to run segmentation of all organelles and get all
    # distance transforms needed for post_processing

    if feature_list is None:
        feature_list = file_check(save_path)

    features = " ".join(feature_list)

    if "cell" in features:
        print("segmenting cells")
        cell_model = org_models.OrganelleModel("cell")
        cell_mask = cell_model.segment(
            img_path=image_path,
            channel=channels["actin"],
            save=True,
            save_path=save_path,
        )

    if "mito" in features:
        print("segmenting mitos")
        mito_model = org_models.OrganelleModel("mito")
        mito_mask = mito_model.segment(
            img_path=image_path,
            channel=channels["mito"],
            save=True,
            save_path=save_path,

        )

    if "lipid" in features:
        print("segmenting lipid droplets")
        lipid_model = org_models.OrganelleModel("lipid_droplet")
        lipid_mask = lipid_model.segment(

            img_path=image_path,
            channel=channels["lipid"],
            save=True,
            save_path=save_path,
        )

    if "perox" in features:
        print("segmenting peroxisomes")
        peroxisome_model = org_models.OrganelleModel("peroxisome")
        peroxisome_mask = peroxisome_model.segment(
            img_path=image_path,
            channel=channels["peroxi"],
            save=True,
            save_path=save_path,
        )

    if "nuclei" in features:
        print("segmenting nuclei")
        nuclei_model = org_models.OrganelleModel("nuclei")
        nuclei_mask = nuclei_model.segment(
            img_path=image_path,
            channel=channels["nuclei"],
            save=True,
            save_path=save_path,
        )

    if "central" in features or "portal" in features:
        print("calculating distance to central and portal veins")
        vein_distance = vein_dist(f'{save_path}/cell_mask.npy', save_path)

    if "bound" in features:
        print("calculating distance to cell boundary")
        cell_edge_distance(save_path)

    return


def file_check(path):
    # check if any files are missing

    available_files = glob.glob(f"{path}/*")

    # Need to add mask files for any new organelles
    required_np_files = [
        "cell_mask.npy",
        "mito_mask.npy",
        "lipid_droplet_mask.npy",
        "peroxisome_mask.npy",
        "nuclei_mask.npy",
        "central_dist.npy",
        "portal_dist.npy",
        "boundary_dist.npy",
    ]

    missing = False
    to_process = []
    for file in required_np_files:
        file_path = f"{path}/" + file
        if file_path not in available_files:
            missing = True
            print(f"missing {file}")
            to_process.append(file)

    if missing is False:
        print("found all required files :)")
    print("-" * 30)

    return to_process


def vein_distance(save_path, vein):
    pv_image = imread(f"{save_path}/{vein}v_distance.tif")
    pv_distance = distance_transform_edt(pv_image)
    np.save(f"{save_path}/{vein}v_distance.npy", pv_distance)

    return pv_distance


def cell_edge_distance(save_path):
    cell_mask = np.load(f"{save_path}/cell_mask.npy")
    outlines = utils.masks_to_outlines(cell_mask)
    dist_transform = distance_transform_edt((outlines == False) * 1)

    np.save(f"{save_path}/boundary_dist.npy", dist_transform)
    return


if __name__ == "__main__":
    test_img_path = "../test_set/0-Actin_Dendra_LD-1_Merge_WD_M1_12pm_dendra2_actin647_lipitox594.lif - TileScan 2_Merged-2.tif"
    path = "../test_set"
    a = file_check(path)
    preprocess(test_img_path, path, a)
