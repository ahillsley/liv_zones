import pandas as pd
import numpy as np

from skimage.measure import regionprops_table

from liv_zones import organelle as org


mito_properties = [
    "mito_density",
    "mito_avg_area",
    "mito_percent_total_area",
    "mito_distance_from_edge",
    "mito_aspect_ratio",
    "mito_solidity",
    "type_1_mito_density",
    "type_2_mito_density",
    "type_3_mito_density",
    "type_1_mito_avg_area",
    "type_2_mito_avg_area",
    "type_3_mito_avg_area",
    "type_1_mito_avg_aspect_ratio",
    "type_2_mito_avg_aspect_ratio",
    "type_3_mito_avg_aspect_ratio",
    "type_1_mito_percent_total_area",
    "type_2_mito_percent_total_area",
    "type_3_mito_percent_total_area",
    "percent_type_1_mito",
    "percent_type_2_mito",
    "percent_type_3_mito",
    "type_1_mito_dist_from_edge",
    "type_2_mito_dist_from_edge",
    "type_3_mito_dist_from_edge",
]

lipid_droplet_properties = [
    "ld_density",
    "ld_avg_area",
    "ld_percent_total_area",
    "ld_distance_from_edge",
    "type_1_ld_density",
    "type_2_ld_density",
    "type_3_ld_density",
    "type_1_ld_avg_area",
    "type_2_ld_avg_area",
    "type_3_ld_avg_area",
    "type_1_ld_percent_total_area",
    "type_2_ld_percent_total_area",
    "type_3_ld_percent_total_area",
    "percent_type_1_ld",
    "percent_type_2_ld",
    "percent_type_3_ld",
    "type_1_ld_dist_from_edge",
    "type_2_ld_dist_from_edge",
    "type_3_ld_dist_from_edge",
]

peroxisome_properties = [
    "peroxisome_density",
    "peroxisome_avg_area",
    "peroxisome_percent_total_area",
    "peroxisome_distance_from_edge",
    "peroxisome_aspect_ratio",
    "peroxisome_solidity",
    "type_1_peroxisome_density",
    "type_2_peroxisome_density",
    "type_3_peroxisome_density",
    "type_1_peroxisome_avg_area",
    "type_2_peroxisome_avg_area",
    "type_3_peroxisome_avg_area",
    "type_1_peroxisome_avg_aspect_ratio",
    "type_2_peroxisome_avg_aspect_ratio",
    "type_3_peroxisome_avg_aspect_ratio",
    "type_1_peroxisome_percent_total_area",
    "type_2_peroxisome_percent_total_area",
    "type_3_peroxisome_percent_total_area",
    "percent_type_1_peroxisome",
    "percent_type_2_peroxisome",
    "percent_type_3_peroxisome",
    "type_1_peroxisome_dist_from_edge",
    "type_2_peroxisome_dist_from_edge",
    "type_3_peroxisome_dist_from_edge",
]

nuclei_properties = [
    "nuclei_density",
    "nuclei_avg_area",
    "nuclei_percent_total_area",
    "nuclei_distance_from_edge",
    "nuclei_aspect_ratio",
]


def cell_features(path, scale):

    property_list = []

    try:
        mito_props = pd.read_csv(f"{path}/mitochondria_properties.csv")
        property_list += mito_properties
    except:
        mito_props = None

    try:
        lipid_props = pd.read_csv(f"{path}/lipid_droplet_properties.csv")
        property_list += lipid_droplet_properties
    except:
        lipid_props = None

    try:
        peroxisome_props = pd.read_csv(f"{path}/peroxisome_properties.csv")
        property_list += peroxisome_properties
    except:
        peroxisome_props = None

    try:
        nuclei_props = pd.read_csv(f"{path}/nuclei_properties.csv")
        property_list += nuclei_properties
    except:
        nuclei_props = None

    masks = org.Masks(path)

    cell_props(
        mito_props=mito_props,
        lipid_props=lipid_props,
        peroxisome_props=peroxisome_props,
        nuclei_props=nuclei_props,
        save_path=path,
        masks=masks,
        scale=scale,
        properties_list=property_list,
    )


def cell_props(
    mito_props,
    lipid_props,
    peroxisome_props,
    nuclei_props,
    save_path,
    masks,
    scale,
    properties_list,
    save=True,
):

    cell_props_dict = regionprops_table(
        masks.cell_mask, properties=("area", "centroid", "label")
    )

    cell_props_1 = pd.DataFrame.from_dict(cell_props_dict)

    cell_props_2 = pd.DataFrame(
        columns=properties_list, index=range(0, len(cell_props_dict["label"]))
    )

    cell_props = pd.concat((cell_props_1, cell_props_2), axis=1)
    cell_props["area"] = cell_props["area"] / scale ** 2
    cell_props["ascini_position"] = org.ascini_position(
        cell_props, masks.cv_distance, masks.pv_distance
    )

    organelles = OrganelleFuncs(index=1, cell_props=cell_props)

    for cell in cell_props["label"]:
        index = cell - 1
        organelles.index = index

        # Mitochondria Properties
        # -----------------------
        if mito_props is not None:
            organelles.single_cell_mitos = mito_props[mito_props["cell_id"] == cell]

            cell_props["mito_density"][index] = organelles.density("mito")
            cell_props["mito_avg_area"][index] = organelles.avg_area("mito")
            cell_props["mito_percent_total_area"][
                index
            ] = organelles.percent_total_area("mito")
            cell_props["mito_aspect_ratio"][index] = organelles.aspect_ratio("mito")
            cell_props["mito_solidity"][index] = organelles.solidity("mito")
            cell_props["mito_distance_from_edge"][
                index
            ] = organelles.distance_from_edge("mito")

            cell_props["type_1_mito_density"][index] = organelles.type_density(
                "mito", 1
            )
            cell_props["type_2_mito_density"][index] = organelles.type_density(
                "mito", 2
            )
            cell_props["type_3_mito_density"][index] = organelles.type_density(
                "mito", 3
            )

            cell_props["type_1_mito_avg_area"][index] = organelles.type_avg_area(
                "mito", 1
            )
            cell_props["type_2_mito_avg_area"][index] = organelles.type_avg_area(
                "mito", 2
            )
            cell_props["type_3_mito_avg_area"][index] = organelles.type_avg_area(
                "mito", 3
            )

            cell_props["type_1_mito_avg_aspect_ratio"][
                index
            ] = organelles.type_aspect_ratio("mito", 1)
            cell_props["type_2_mito_avg_aspect_ratio"][
                index
            ] = organelles.type_aspect_ratio("mito", 2)
            cell_props["type_3_mito_avg_aspect_ratio"][
                index
            ] = organelles.type_aspect_ratio("mito", 3)

            cell_props["type_1_mito_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("mito", 1)
            cell_props["type_2_mito_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("mito", 2)
            cell_props["type_3_mito_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("mito", 3)

            cell_props["percent_type_1_mito"][
                index
            ] = organelles.type_percent_organelles("mito", 1)
            cell_props["percent_type_2_mito"][
                index
            ] = organelles.type_percent_organelles("mito", 2)
            cell_props["percent_type_3_mito"][
                index
            ] = organelles.type_percent_organelles("mito", 3)

            cell_props["type_1_mito_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("mito", 1)
            cell_props["type_2_mito_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("mito", 2)
            cell_props["type_3_mito_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("mito", 3)

        # Lipid Droplet Properties
        # ------------------------
        if lipid_props is not None:
            organelles.single_cell_lipid_droplets = lipid_props[
                lipid_props["cell_id"] == cell
            ]

            cell_props["ld_density"][index] = organelles.density("ld")
            cell_props["ld_avg_area"][index] = organelles.avg_area("ld")
            cell_props["ld_percent_total_area"][index] = organelles.percent_total_area(
                "ld"
            )
            cell_props["ld_distance_from_edge"][index] = organelles.distance_from_edge(
                "ld"
            )

            cell_props["type_1_ld_density"][index] = organelles.type_density("ld", 1)
            cell_props["type_2_ld_density"][index] = organelles.type_density("ld", 2)
            cell_props["type_3_ld_density"][index] = organelles.type_density("ld", 3)

            cell_props["type_1_ld_avg_area"][index] = organelles.type_avg_area("ld", 1)
            cell_props["type_2_ld_avg_area"][index] = organelles.type_avg_area("ld", 2)
            cell_props["type_3_ld_avg_area"][index] = organelles.type_avg_area("ld", 3)

            cell_props["type_1_ld_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("ld", 1)
            cell_props["type_2_ld_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("ld", 2)
            cell_props["type_3_ld_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("ld", 3)

            cell_props["percent_type_1_ld"][index] = organelles.type_percent_organelles(
                "ld", 1
            )
            cell_props["percent_type_2_ld"][index] = organelles.type_percent_organelles(
                "ld", 2
            )
            cell_props["percent_type_3_ld"][index] = organelles.type_percent_organelles(
                "ld", 3
            )

            cell_props["type_1_ld_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("ld", 1)
            cell_props["type_2_ld_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("ld", 2)
            cell_props["type_3_ld_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("ld", 3)

        # Peroxisome Properties
        # ---------------------
        if peroxisome_props is not None:
            organelles.single_cell_peroxisomes = peroxisome_props[
                peroxisome_props["cell_id"] == cell
            ]

            cell_props["peroxisome_density"][index] = organelles.density("peroxi")
            cell_props["peroxisome_avg_area"][index] = organelles.avg_area("peroxi")
            cell_props["peroxisome_percent_total_area"][
                index
            ] = organelles.percent_total_area("peroxi")
            cell_props["peroxisome_aspect_ratio"][index] = organelles.aspect_ratio(
                "peroxi"
            )
            cell_props["peroxisome_solidity"][index] = organelles.solidity("peroxi")
            cell_props["peroxisome_distance_from_edge"][
                index
            ] = organelles.distance_from_edge("peroxi")

            cell_props["type_1_peroxisome_density"][index] = organelles.type_density(
                "peroxi", 1
            )
            cell_props["type_2_peroxisome_density"][index] = organelles.type_density(
                "peroxi", 2
            )
            cell_props["type_3_peroxisome_density"][index] = organelles.type_density(
                "peroxi", 3
            )

            cell_props["type_1_peroxisome_avg_area"][index] = organelles.type_avg_area(
                "peroxi", 1
            )
            cell_props["type_2_peroxisome_avg_area"][index] = organelles.type_avg_area(
                "peroxi", 2
            )
            cell_props["type_3_peroxisome_avg_area"][index] = organelles.type_avg_area(
                "peroxi", 3
            )

            cell_props["type_1_peroxisome_avg_aspect_ratio"][
                index
            ] = organelles.type_aspect_ratio("peroxi", 1)
            cell_props["type_2_peroxisome_avg_aspect_ratio"][
                index
            ] = organelles.type_aspect_ratio("peroxi", 2)
            cell_props["type_3_peroxisome_avg_aspect_ratio"][
                index
            ] = organelles.type_aspect_ratio("peroxi", 3)

            cell_props["type_1_peroxisome_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("peroxi", 1)
            cell_props["type_2_peroxisome_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("peroxi", 2)
            cell_props["type_3_peroxisome_percent_total_area"][
                index
            ] = organelles.type_percent_total_area("peroxi", 3)

            cell_props["percent_type_1_peroxisome"][
                index
            ] = organelles.type_percent_organelles("peroxi", 1)
            cell_props["percent_type_2_peroxisome"][
                index
            ] = organelles.type_percent_organelles("peroxi", 2)
            cell_props["percent_type_3_peroxisome"][
                index
            ] = organelles.type_percent_organelles("peroxi", 3)

            cell_props["type_1_peroxisome_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("peroxi", 1)
            cell_props["type_2_peroxisome_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("peroxi", 2)
            cell_props["type_3_peroxisome_dist_from_edge"][
                index
            ] = organelles.type_dist_from_edge("peroxi", 3)

        if nuclei_props is not None:
            organelles.single_cell_nuclei = nuclei_props[
                nuclei_props["cell_id"] == cell
            ]

            cell_props["nuclei_density"][index] = organelles.density("nuclei")
            cell_props["nuclei_avg_area"][index] = organelles.avg_area("nuclei")
            cell_props["nuclei_percent_total_area"][
                index
            ] = organelles.percent_total_area("nuclei")
            cell_props["nuclei_aspect_ratio"][index] = organelles.aspect_ratio("nuclei")
            cell_props["nuclei_distance_from_edge"][
                index
            ] = organelles.distance_from_edge("nuclei")

    if save is True:
        cell_props.to_csv(f"{save_path}/average_properties_per_cell.csv")

    return cell_props


class OrganelleFuncs:
    """
    Functions to average the properties of individual organelles at the
    single-cell level

    Functions take up to 2 arguments:

        organelle (str):
            - valid options are 'mito' (mitochondria) or 'ld' (lipid_droplet)
            - specifies which organelle to apply function to

        org_type (int):
            - valid options are 1, 2, 3
            - specifies which organelle sub-classification to consider
            - only valid for select funtions, depending on which properties
                have been sub-classified

    """

    def __init__(self, index, cell_props):

        self.index = index
        self.cell_props = cell_props
        self.single_cell_mitos = None
        self.single_cell_lipid_droplets = None
        self.single_cell_peroxisomes = None
        self.single_cell_nuclei = None

    def density(self, organelle, org_type=False):
        if organelle == "mito":
            organelle_count = len(self.single_cell_mitos)

        elif organelle == "ld":
            organelle_count = len(self.single_cell_lipid_droplets)

        elif organelle == "peroxi":
            organelle_count = len(self.single_cell_peroxisomes)

        elif organelle == "nuclei":
            organelle_count = len(self.single_cell_nuclei)

        return organelle_count / self.cell_props["area"][self.index]

    def avg_area(self, organelle):
        if organelle == "mito":
            return np.mean(self.single_cell_mitos["area"])

        elif organelle == "ld":
            return np.mean(self.single_cell_lipid_droplets["area"])

        elif organelle == "peroxi":
            return np.mean(self.single_cell_peroxisomes["area"])

        elif organelle == "nuclei":
            return np.mean(self.single_cell_nuclei["area"])

    def percent_total_area(self, organelle):
        if organelle == "mito":
            organelle_area = np.sum(self.single_cell_mitos["area"])

        elif organelle == "ld":
            organelle_area = np.sum(self.single_cell_lipid_droplets["area"])

        elif organelle == "peroxi":
            organelle_area = np.sum(self.single_cell_peroxisomes["area"])

        elif organelle == "nuclei":
            organelle_area = np.sum(self.single_cell_nuclei["area"])

        return organelle_area / self.cell_props["area"][self.index]

    def aspect_ratio(self, organelle):
        if organelle == "mito":
            return np.mean(self.single_cell_mitos["aspect_ratio"])

        elif organelle == "peroxi":
            return np.mean(self.single_cell_peroxisomes["aspect_ratio"])

        elif organelle == "nuclei":
            return np.mean(self.single_cell_nuclei["aspect_ratio"])

    def solidity(self, organelle):
        if organelle == "mito":
            return np.mean(self.single_cell_mitos["solidity"])

        elif organelle == "peroxi":
            return np.mean(self.single_cell_peroxisomes["solidity"])

    def distance_from_edge(self, organelle):
        if organelle == "mito":
            organelle_type = self.single_cell_mitos

        elif organelle == "ld":
            organelle_type = self.single_cell_lipid_droplets

        elif organelle == "peroxi":
            organelle_type = self.single_cell_peroxisomes

        elif organelle == "nuclei":
            organelle_type = self.single_cell_nuclei

        return dist_from_edge(self.cell_props.iloc[self.index], organelle_type)

    def type_density(self, organelle, org_type):
        if organelle == "mito":
            organelle_count = np.sum(self.single_cell_mitos[f"aspect_type_{org_type}"])

        elif organelle == "ld":
            organelle_count = np.sum(
                self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            )

        elif organelle == "peroxi":
            organelle_count = np.sum(
                self.single_cell_peroxisomes[f"aspect_type_{org_type}"]
            )

        return organelle_count / self.cell_props["area"][self.index]

    def type_avg_area(self, organelle, org_type):
        if organelle == "mito":
            subset = self.single_cell_mitos[f"aspect_type_{org_type}"]
            return np.mean(self.single_cell_mitos["area"][subset])

        elif organelle == "ld":
            subset = self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            return np.mean(self.single_cell_lipid_droplets["area"][subset])

        elif organelle == "peroxi":
            subset = self.single_cell_peroxisomes[f"aspect_type_{org_type}"]
            return np.mean(self.single_cell_peroxisomes["area"][subset])

    def type_aspect_ratio(self, organelle, org_type):
        if organelle == "mito":
            subset = self.single_cell_mitos[f"aspect_type_{org_type}"]
            return np.mean(self.single_cell_mitos["aspect_ratio"][subset])

        if organelle == "peroxi":
            subset = self.single_cell_mitos[f"aspect_type_{org_type}"]
            return np.mean(self.single_cell_mitos["aspect_ratio"][subset])

    def type_percent_total_area(self, organelle, org_type):
        if organelle == "mito":
            organelle_area = np.sum(
                self.single_cell_mitos["area"][
                    self.single_cell_mitos[f"aspect_type_{org_type}"]
                ]
            )

        elif organelle == "ld":
            organelle_area = np.sum(
                self.single_cell_lipid_droplets["area"][
                    self.single_cell_lipid_droplets[f"area_type_{org_type}"]
                ]
            )

        elif organelle == "peroxi":
            organelle_area = np.sum(
                self.single_cell_peroxisomes["area"][
                    self.single_cell_peroxisomes[f"aspect_type_{org_type}"]
                ]
            )

        return organelle_area / self.cell_props["area"][self.index]

    def type_percent_organelles(self, organelle, org_type):
        if organelle == "mito":
            return np.sum(self.single_cell_mitos[f"aspect_type_{org_type}"]) / len(
                self.single_cell_mitos
            )

        elif organelle == "ld":
            return np.sum(
                self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            ) / len(self.single_cell_lipid_droplets)

        elif organelle == "peroxi":
            return np.sum(self.single_cell_peroxisomes[f"aspect_type_{org_type}"]) / len(
                self.single_cell_peroxisomes
            )

    def type_dist_from_edge(self, organelle, org_type):
        if organelle == "mito":
            subset = self.single_cell_mitos[
                self.single_cell_mitos[f"aspect_type_{org_type}"]
            ]

        elif organelle == "ld":
            subset = self.single_cell_lipid_droplets[
                self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            ]

        elif organelle == "peroxi":
            subset = self.single_cell_peroxisomes[
                self.single_cell_peroxisomes[f"aspect_type_{org_type}"]
            ]

        return dist_from_edge(self.cell_props.iloc[self.index], subset)


def dist_from_edge(cell, organelle_list):
    cc_x = float(cell["centroid-0"])
    cc_y = float(cell["centroid-1"])

    o_s = np.asarray([organelle_list["centroid-0"], organelle_list["centroid-1"]])

    cc_s = np.asarray([np.repeat(cc_x, o_s.shape[1]), np.repeat(cc_y, o_s.shape[1])])
    # in pixels
    dist_to_center = np.linalg.norm(cc_s - o_s, axis=0)
    dist_to_edge = np.asarray(organelle_list["boundry_dist"])
    rel_position = dist_to_center / (dist_to_center + dist_to_edge)

    return np.mean(rel_position)
