import numpy as np
import pandas as pd

from skimage.measure import regionprops_table


# Adjust default values here
# --------------------------
mito_aspect_split = (1.352, 1.970)
ld_area_split = (0.88, 2.283)
peroxisome_aspect_split = (1.327, 1.843)  # completely made up, need to change


class Masks:
    def __init__(self, path):
        self.cell_mask = np.load(f"{path}/cell_mask.npy")
        self.cv_distance = np.load(f"{path}/central_dist.npy")
        self.pv_distance = np.load(f"{path}/portal_dist.npy")
        self.cell_edge_distance = np.load(f"{path}/boundary_dist.npy")


def organelle_features(
    path,
    scale,
    mito_aspect_split=mito_aspect_split,
    ld_area_split=ld_area_split,
    organelle_list=["mitos", "lipid_droplets", "peroxisomes", "nuclei"],
):

    cell_level_masks = Masks(path)

    organelles = " ".join(organelle_list)

    if "mito" in organelles:
        mito_mask = np.load(f"{path}/mito_mask.npy")
        mito_properties(mito_mask, cell_level_masks, path, scale, mito_aspect_split)

    if "lipid" in organelles:
        lipid_droplet_mask = np.load(f"{path}/lipid_droplet_mask.npy")
        lipid_droplet_properties(
            lipid_droplet_mask, cell_level_masks, path, scale, ld_area_split
        )

    if "perox" in organelles:
        peroxisome_mask = np.load(f"{path}/peroxisome_mask.npy")
        peroxisome_properties(
            peroxisome_mask, cell_level_masks, path, scale, peroxisome_aspect_split
        )

    if "nuclei" in organelles:
        nuclei_mask = np.load(f"{path}/nuclei_mask.npy")
        nuclei_properties(
                nuclei_mask, cell_level_masks, path, scale
        )

    return


def mito_properties(
    mito_mask, masks, save_path, scale, mito_aspect_split, save=True  # pixels / micron
):

    # returns props in pixel units
    mito_props_dict = regionprops_table(
        mito_mask,
        properties=(
            "label",
            "area",
            "perimeter",
            "centroid",
            "axis_major_length",
            "axis_minor_length",
            "solidity",
        ),
    )

    mito_props = pd.DataFrame.from_dict(mito_props_dict)

    # convert from pixels to microns
    mito_props["area"] = mito_props["area"] / (scale ** 2)
    mito_props["perimeter"] = mito_props["perimeter"] / scale
    mito_props["axis_major_length"] = mito_props["axis_major_length"] / scale
    mito_props["axis_minor_length"] = mito_props["axis_minor_length"] / scale

    mito_props["cell_id"] = map_to_cell(mito_props, masks.cell_mask)
    mito_props["aspect_ratio"] = (
        mito_props["axis_major_length"] / mito_props["axis_minor_length"]
    )

    # use formula to calc boundry to cell edge
    mito_props["boundry_dist"] = map_to_cell(mito_props, masks.cell_edge_distance)

    # use new ascini distance measure
    # portal vein = 1, central vein = -1
    mito_props["ascini_position"] = ascini_position(
        mito_props, masks.cv_distance, masks.pv_distance
    )

    (
        mito_props["aspect_type_1"],
        mito_props["aspect_type_2"],
        mito_props["aspect_type_3"],
    ) = split_types(mito_props, "aspect_ratio", mito_aspect_split)
    if save is True:
        mito_props.to_csv(f"{save_path}/mitochondria_properties.csv")

    return mito_props


def lipid_droplet_properties(
    lipid_droplet_mask,
    masks,
    save_path,
    scale,  # pixels per micron
    ld_area_split,
    save=True,
):

    ld_props_dict = regionprops_table(
        lipid_droplet_mask,
        properties=(
            "label",
            "area",
            "perimeter",
            "centroid",
            "axis_major_length",
            "axis_minor_length",
        ),
    )

    ld_props = pd.DataFrame.from_dict(ld_props_dict)

    # convert from pixels to microns
    ld_props["area"] = ld_props["area"] / (scale ** 2)
    ld_props["perimeter"] = ld_props["perimeter"] / scale
    ld_props["axis_major_length"] = ld_props["axis_major_length"] / scale
    ld_props["axis_minor_length"] = ld_props["axis_minor_length"] / scale

    ld_props["cell_id"] = map_to_cell(ld_props, masks.cell_mask)
    ld_props["aspect_ratio"] = (
        ld_props["axis_major_length"] / ld_props["axis_minor_length"]
    )

    ld_props["boundry_dist"] = map_to_cell(ld_props, masks.cell_edge_distance)

    ld_props["ascini_position"] = ascini_position(
        ld_props, masks.cv_distance, masks.pv_distance
    )

    (
        ld_props["area_type_1"],
        ld_props["area_type_2"],
        ld_props["area_type_3"],
    ) = split_types(ld_props, "area", ld_area_split)

    if save is True:
        ld_props.to_csv(f"{save_path}/lipid_droplet_properties.csv")

    return ld_props


def peroxisome_properties(
    peroxisome_mask,
    masks,
    save_path,
    scale,
    peroxisome_aspect_split,
    save=True,  # pixels / micron
):

    # returns props in pixel units
    peroxi_props_dict = regionprops_table(
        peroxisome_mask,
        properties=(
            "label",
            "area",
            "perimeter",
            "centroid",
            "axis_major_length",
            "axis_minor_length",
            "solidity",
        ),
    )

    peroxi_props = pd.DataFrame.from_dict(peroxi_props_dict)

    # convert from pixels to microns
    peroxi_props["area"] = peroxi_props["area"] / (scale ** 2)
    peroxi_props["perimeter"] = peroxi_props["perimeter"] / scale
    peroxi_props["axis_major_length"] = peroxi_props["axis_major_length"] / scale
    peroxi_props["axis_minor_length"] = peroxi_props["axis_minor_length"] / scale

    peroxi_props["cell_id"] = map_to_cell(peroxi_props, masks.cell_mask)
    peroxi_props["aspect_ratio"] = (
        peroxi_props["axis_major_length"] / peroxi_props["axis_minor_length"]
    )

    # use formula to calc boundry to cell edge
    peroxi_props["boundry_dist"] = map_to_cell(peroxi_props, masks.cell_edge_distance)

    # use new ascini distance measure
    # portal vein = 1, central vein = -1
    peroxi_props["ascini_position"] = ascini_position(
        peroxi_props, masks.cv_distance, masks.pv_distance
    )

    (
        peroxi_props["aspect_type_1"],
        peroxi_props["aspect_type_2"],
        peroxi_props["aspect_type_3"],
    ) = split_types(peroxi_props, "aspect_ratio", peroxisome_aspect_split)
    if save is True:
        peroxi_props.to_csv(f"{save_path}/peroxisome_properties.csv")

    return peroxi_props


def nuclei_properties(nuclei_mask, masks, save_path, scale, save=True):

    # returns props in pixel units
    nuclei_props_dict = regionprops_table(
        nuclei_mask,
        properties=(
            "label",
            "area",
            "perimeter",
            "centroid",
            "axis_major_length",
            "axis_minor_length",
        ),
    )

    nuclei_props = pd.DataFrame.from_dict(nuclei_props_dict)

    # convert from pixels to microns
    nuclei_props["area"] = nuclei_props["area"] / (scale ** 2)
    nuclei_props["perimeter"] = nuclei_props["perimeter"] / scale
    nuclei_props["axis_major_length"] = nuclei_props["axis_major_length"] / scale
    nuclei_props["axis_minor_length"] = nuclei_props["axis_minor_length"] / scale

    nuclei_props["cell_id"] = map_to_cell(nuclei_props, masks.cell_mask)
    nuclei_props["aspect_ratio"] = (
        nuclei_props["axis_major_length"] / nuclei_props["axis_minor_length"]
    )

    # use formula to calc boundry to cell edge
    nuclei_props["boundry_dist"] = map_to_cell(nuclei_props, masks.cell_edge_distance)

    # use new ascini distance measure
    # portal vein = 1, central vein = -1
    nuclei_props["ascini_position"] = ascini_position(
        nuclei_props, masks.cv_distance, masks.pv_distance
    )

    if save is True:
        nuclei_props.to_csv(f"{save_path}/nuclei_properties.csv")

    return nuclei_props


def map_to_cell(props, mask):
    y_cords = np.asarray(props["centroid-0"], dtype="int")
    x_cords = np.asarray(props["centroid-1"], dtype="int")

    return mask[y_cords, x_cords]


def split_types(organelle_list, prop, bounds):
    """
    Split the given dataframe into types based on "prop" and cutoffs
    """

    type_1 = organelle_list[prop] < bounds[0]
    type_3 = organelle_list[prop] >= bounds[1]
    type_2 = np.logical_or(type_1, type_3) == False

    return type_1, type_2, type_3


def ascini_position(props, cv_distance, pv_distance):
    cv_distance = map_to_cell(props, cv_distance)
    pv_distance = map_to_cell(props, pv_distance)
    ascini_position = (pv_distance - cv_distance) / (pv_distance + cv_distance)

    return ascini_position
