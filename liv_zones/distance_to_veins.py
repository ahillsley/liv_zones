# %%
from skimage.morphology import remove_small_objects
from skimage.measure import label, regionprops
import numpy as np
import matplotlib.pyplot as plt
from skimage.draw import ellipse
from scipy import optimize
from scipy.ndimage import distance_transform_edt
#import typer


def get_veins(data):
    holes = (data == 0).astype(int)
    labelled_holes = label(holes)
    # Get sizes of objects in labelled image
    object_sizes = sorted(np.bincount(labelled_holes.ravel()), reverse=True)
    # Remove all object smaller than the second largest
    vein_size = object_sizes[2]
    veins = remove_small_objects(labelled_holes, vein_size)
    return veins


def cost(params, mask):
    """
    Cost is the number of pixels in the ellipse that are not in the object + the number of pixels in the object that are not in the ellipse

    mask: binary image of object
    x0, y0: center of ellipse
    a, b: major and minor axis lengths
    theta: angle of rotation
    """
    x0, y0, a, b, theta = params
    ellipse_fit = ellipse(x0, y0, a, b, rotation=theta, shape=mask.shape)
    ellipse_mask = np.zeros_like(mask)
    ellipse_mask[ellipse_fit] = 1
    # Minimize the number of 0 values within the ellipse
    # Minimize the number of 1 values outside the ellipse
    return np.sum((ellipse_mask == 0) & (mask == 1)) + np.sum(
        (ellipse_mask == 1) & (mask == 0)
    )


def fit_ellipses(data):
    """
    For every object in the data, fit an ellipse to it such that it maximises true parts of the object while minimizing background inside the ellipse.
    """
    props = regionprops(data)

    fits = {}
    for obj in props:
        # Make a binary mask of the object
        mask = (data == obj.label).astype(int)
        x0, y0 = obj.centroid
        a, b = obj.axis_major_length, obj.axis_minor_length
        theta = 0  # obj.orientation
        # Optimize the cost function
        x0, y0, a, b, theta = optimize.fmin(
            cost, x0=(x0, y0, a, b, theta), args=(mask,), disp=True
        )
        # Create a final fit that has the same shape as the original image but is masked by the ellipse
        ellipse_fit = ellipse(x0, y0, a, b, rotation=theta, shape=data.shape)
        fits[y0] = ellipse_fit
    # Rename fits such that the lower value of x0 is called "central" and the other is "portal"
    central, portal = sorted(fits.keys())
    fits["central"] = fits.pop(central)
    fits["portal"] = fits.pop(portal)
    return fits


def main(input_file, output_dir):
    """
    Generates distance transform images for the central and portal vein.

    Args:
        input_file: path to the cell masks
        output_dir: path to the output directory
            The output are saved as "central_dist.npy" and "portal_dist.npy" in the output directory.

    Returns:
        None
    """
    data = np.load(input_file)
    veins = get_veins(data)
    fits = fit_ellipses(veins)
    for name, fit in fits.items():
        img = np.ones_like(veins)
        img[fit] = 0
        dist = distance_transform_edt(img)
        # Save
        np.save(f"{output_dir}/{name}_dist.npy", dist)


if __name__ == "__main__":
    #typer.run(main)
    pass

