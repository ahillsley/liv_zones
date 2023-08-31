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

property_list = mito_properties + lipid_droplet_properties


def cell_features(path, scale):

    lipid_props = pd.read_csv(f"{path}/lipid_droplet_properties.csv")
    mito_props = pd.read_csv(f"{path}/mitochondria_properties.csv")

    masks = org.Masks(path)

    cell_props(
        mito_props=mito_props,
        lipid_props=lipid_props,
        save_path=path,
        masks=masks,
        scale=scale,
    )


def cell_props(
    mito_props,
    lipid_props,
    save_path,
    masks,
    scale,
    properties_list=property_list,
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
    cell_props["area"] = cell_props["area"] / scale**2
    cell_props["ascini_position"] = org.ascini_position(
        cell_props, masks.cv_distance, masks.pv_distance
    )

    organelles = OrganelleFuncs(
        index=1, mito_props=mito_props, lipid_props=lipid_props, cell_props=cell_props
    )

    for cell in cell_props["label"]:
        index = cell - 1
        organelles.index = index
        organelles.single_cell_mitos = mito_props[mito_props["cell_id"] == cell]
        organelles.single_cell_lipid_droplets = lipid_props[
            lipid_props["cell_id"] == cell
        ]

        # Mitochondria Properties
        # -----------------------
        cell_props["mito_density"][index] = organelles.density("mito")
        cell_props["mito_avg_area"][index] = organelles.avg_area("mito")
        cell_props["mito_percent_total_area"][index] = organelles.percent_total_area(
            "mito"
        )
        cell_props["mito_aspect_ratio"][index] = organelles.mito_aspect_ratio()
        cell_props["mito_solidity"][index] = organelles.mito_solidity()
        cell_props["mito_distance_from_edge"][index] = organelles.distance_from_edge(
            "mito"
        )

        cell_props["type_1_mito_density"][index] = organelles.type_density("mito", 1)
        cell_props["type_2_mito_density"][index] = organelles.type_density("mito", 2)
        cell_props["type_3_mito_density"][index] = organelles.type_density("mito", 3)

        cell_props["type_1_mito_avg_area"][index] = organelles.type_avg_area("mito", 1)
        cell_props["type_2_mito_avg_area"][index] = organelles.type_avg_area("mito", 2)
        cell_props["type_3_mito_avg_area"][index] = organelles.type_avg_area("mito", 3)

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

        cell_props["percent_type_1_mito"][index] = organelles.type_percent_organelles(
            "mito", 1
        )
        cell_props["percent_type_2_mito"][index] = organelles.type_percent_organelles(
            "mito", 2
        )
        cell_props["percent_type_3_mito"][index] = organelles.type_percent_organelles(
            "mito", 3
        )

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
        cell_props["ld_density"][index] = organelles.density("ld")
        cell_props["ld_avg_area"][index] = organelles.avg_area("ld")
        cell_props["ld_percent_total_area"][index] = organelles.percent_total_area("ld")
        cell_props["ld_distance_from_edge"][index] = organelles.distance_from_edge("ld")

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

        cell_props["type_1_ld_dist_from_edge"][index] = organelles.type_dist_from_edge(
            "ld", 1
        )
        cell_props["type_2_ld_dist_from_edge"][index] = organelles.type_dist_from_edge(
            "ld", 2
        )
        cell_props["type_3_ld_dist_from_edge"][index] = organelles.type_dist_from_edge(
            "ld", 3
        )

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

    def __init__(self, index, mito_props, lipid_props, cell_props):

        self.index = index
        self.cell_props = cell_props
        self.single_cell_mitos = None
        self.single_cell_lipid_droplets = None

    def density(self, organelle, org_type=False):
        if organelle == "mito":
            organelle_count = len(self.single_cell_mitos)

        elif organelle == "ld":
            organelle_count = len(self.single_cell_lipid_droplets)

        return organelle_count / self.cell_props["area"][self.index]

    def avg_area(self, organelle):
        if organelle == "mito":
            return np.mean(self.single_cell_mitos["area"])

        elif organelle == "ld":
            return np.mean(self.single_cell_lipid_droplets["area"])

    def percent_total_area(self, organelle):
        if organelle == "mito":
            organelle_area = np.sum(self.single_cell_mitos["area"])

        elif organelle == "ld":
            organelle_area = np.sum(self.single_cell_lipid_droplets["area"])

        return organelle_area / self.cell_props["area"][self.index]

    def mito_aspect_ratio(self):
        return np.mean(self.single_cell_mitos["aspect_ratio"])

    def mito_solidity(self):
        return np.mean(self.single_cell_mitos["solidity"])

    def distance_from_edge(self, organelle):
        if organelle == "mito":
            organelle_type = self.single_cell_mitos

        elif organelle == "ld":
            organelle_type = self.single_cell_lipid_droplets

        return dist_from_edge(self.cell_props.iloc[self.index], organelle_type)

    def type_density(self, organelle, org_type):
        if organelle == "mito":
            organelle_count = np.sum(self.single_cell_mitos[f"aspect_type_{org_type}"])
        elif organelle == "ld":
            organelle_count = np.sum(
                self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            )

        return organelle_count / self.cell_props["area"][self.index]

    def type_avg_area(self, organelle, org_type):
        if organelle == "mito":
            subset = self.single_cell_mitos[f"aspect_type_{org_type}"]
            return np.mean(self.single_cell_mitos["area"][subset])

        elif organelle == "ld":
            subset = self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            return np.mean(self.single_cell_lipid_droplets["area"][subset])

    def type_aspect_ratio(self, organelle, org_type):
        if organelle == "mito":
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

        return organelle_area / self.cell_props["area"][self.index]

    def type_percent_organelles(self, organelle, org_type):
        if organelle == "mito":
            return np.sum(self.single_cell_mitos[f"aspect_type_{org_type}"]) / len(
                self.single_cell_mitos
            )

        if organelle == "ld":
            return np.sum(
                self.single_cell_lipid_droplets[f"area_type_{org_type}"]
            ) / len(self.single_cell_lipid_droplets)

    def type_dist_from_edge(self, organelle, org_type):
        if organelle == "mito":
            subset = self.single_cell_mitos[
                self.single_cell_mitos[f"aspect_type_{org_type}"]
            ]

        if organelle == "ld":
            subset = self.single_cell_lipid_droplets[
                self.single_cell_lipid_droplets[f"area_type_{org_type}"]
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
