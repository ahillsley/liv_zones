
import os
import glob
import math
import pickle
import collections
import re

import numpy as np
from PIL import Image, ImageFilter
import cv2 as cv
from skimage.filters import rank
from scipy import ndimage
from matplotlib import patches
import tifffile
from skimage.morphology import disk

Image.MAX_IMAGE_PIXELS = 500000000


# ================================================================
# Robust filename parsing helpers
# ================================================================
def extract_z_index(filename):
    """
    Extract z index from filenames like:
    - Region 2_Merged_z002_RAW_ch03.tif
    - Region 7_Merged--Z00--C01.tif
    - ..._z12_...
    - ...--Z003--...
    """
    # Accept "_z###", "-z###", "--Z###--", etc. Case-insensitive.
    m = re.search(r"(?:^|[^A-Za-z0-9])z(\d+)(?:[^A-Za-z0-9]|$)", filename, flags=re.IGNORECASE)
    if m is None:
        raise ValueError(f"Could not find z index in filename: {filename}")
    return int(m.group(1))



# ================================================================
# Helper to get reference size from first image
# ================================================================
def getSampleAsinus(organelle_dir):
    """
    Use the first TIFF image in organelle_dir to determine image size.

    Returns
    -------
    sample_asinus : np.ndarray
        2D array of the first image (full frame).
    crop_info : dict
        Dictionary containing 'height' and 'width' (for compatibility).
    """
    sample_image_path = glob.glob(f"{organelle_dir}/*.tif")
    if not sample_image_path:
        raise FileNotFoundError(f"No TIFF images found in {organelle_dir}")

    sample_asinus_slice = Image.open(sample_image_path[0])
    sample_asinus = np.array(sample_asinus_slice)

    # Assume grayscale 2D images
    H, W = sample_asinus.shape
    crop_info = {
        "height": H,
        "width": W,
    }

    return sample_asinus, crop_info


# ================================================================
# DAPI processing – full frame, no cropping
# ================================================================
def getDAPI(nuclei_dir, nSlices, stackNum, crop_info):
    """
    Load DAPI channel as full-frame Z-stack and return max projection (median-filtered).
    """
    H = crop_info["height"]
    W = crop_info["width"]

    zstack = np.zeros((nSlices, H, W), dtype="uint8")
    print("DAPI zstack shape:", zstack.shape)

    files = glob.glob(f"{nuclei_dir}/*.tif")
    # Robust z sorting (supports names with extra tokens like _RAW_)
    files.sort(key=lambda x: extract_z_index(os.path.basename(x)))

    for i in range(nSlices):
        idx = nSlices * stackNum + i
        asinus_slice = Image.open(files[idx])
        asinus = np.array(asinus_slice, dtype="uint8")  # full frame
        zstack[i] = asinus

    dapi_maxproj = np.max(zstack, axis=0)
    dapi_maxproj_filtered = ndimage.median_filter(dapi_maxproj, size=1)
    return dapi_maxproj_filtered


# ================================================================
# Organelles processing – full frame, channel-specific filters
# ================================================================
def getOrgStacks(organelle_dir, channels, nSlices, stackNum, crop_info):
    """
    Build organelle max projections for all channels, full frame, no cropping.

    Returns
    -------
    asinus_maxproj : np.ndarray
        Shape (len(channels), H, W), dtype=uint16.
    """
    H = crop_info["height"]
    W = crop_info["width"]

    asinus_maxproj = np.zeros((len(channels), H, W), dtype="uint16")

    for i, (key, values) in enumerate(channels.items()):
        print(f"Processing channel {key} ({i+1}/{len(channels)})")

        # All TIFFs for this channel
        org_files = glob.glob(f"{organelle_dir}/*{values}.tif")
        # Robust z sorting (supports names with extra tokens like _RAW_)
        org_files.sort(key=lambda x: extract_z_index(os.path.basename(x)))

        # ============================================================
        # CASE 1: glucagon / insulin / actin / lipid
        # ONLY max Z projection, no filters (as in your original)
        # ============================================================
        if key in ["glucagon", "insulin", "actin", "lipid"]:
            zstack = np.zeros((nSlices, H, W), dtype="uint16")

            for j in range(nSlices):
                idx = nSlices * stackNum + j
                asinus_slice = Image.open(org_files[idx])
                asinusa = np.array(asinus_slice)  # full frame
                zstack[j] = asinusa.astype("uint16")

            asinus_maxproj[i] = np.max(zstack, axis=0)

        # ============================================================
        # CASE 2: everything else (mito, peroxi, etc.)
        # For mito/peroxi: pass RAW uint8 planes into corrections()
        # so it can do: blur/subtract -> median(radius=1) -> MIP
        # ============================================================
        else:
            zstack_u8 = np.zeros((nSlices, H, W), dtype="uint8")

            for j in range(nSlices):
                idx = nSlices * stackNum + j
                asinus_slice = Image.open(org_files[idx])
                asinusa = np.array(asinus_slice)  # full frame

                # Force uint8 (matches your original approach via np.uint8(...))
                zstack_u8[j] = asinusa.astype("uint8")

            corrected_2d = corrections(key, zstack_u8, None)

            # Store into uint16 output stack (safe upcast)
            asinus_maxproj[i] = corrected_2d.astype("uint16")

        print(f"Key processed: {key}")

    return asinus_maxproj


# ================================================================
# Corrections – UPDATED mito/peroxi logic
# ================================================================
# ================================================================
# Separate correction functions (mito vs peroxi)
# ================================================================
def corrections_mito(image_u8_stack):
    """
    Mito correction:
      1) Gaussian blur per z-plane (sigma=150, no blur across Z)
      2) Subtract: raw - blurred (clip 0..255)
      3) ImageJ-like median filter per z-plane with disk radius=2
      4) Max projection at the end
    """
    img_u8 = image_u8_stack.astype(np.uint8)

    blur = ndimage.gaussian_filter(img_u8.astype(np.float32), sigma=(0, 150, 150))

    sub = img_u8.astype(np.float32) - blur
    sub = np.clip(sub, 0, 255).astype(np.uint8)

    fp = disk(2)
    sub_med = np.empty_like(sub, dtype=np.uint8)
    for z in range(sub.shape[0]):
        sub_med[z] = rank.median(sub[z], footprint=fp)

    out = np.max(sub_med, axis=0).astype(np.uint8)
    return out


def corrections_peroxi(image_u8_stack):
    """
    Peroxisome correction:
    Currently identical to mito, but kept separate so you can change later
    without touching mitochondria.
    """
    img_u8 = image_u8_stack.astype(np.uint8)

    blur = ndimage.gaussian_filter(img_u8.astype(np.float32), sigma=(0, 150, 150))

    sub = img_u8.astype(np.float32) - blur
    sub = np.clip(sub, 0, 255).astype(np.uint8)

    fp = disk(2)
    sub_med = np.empty_like(sub, dtype=np.uint8)
    for z in range(sub.shape[0]):
        sub_med[z] = rank.median(sub[z], footprint=fp)

    out = np.max(sub_med, axis=0).astype(np.uint8)
    return out


#================================================================
# Output Writer
# ================================================================
def outputTIFs(root_dir, stackNum, asinus_maxproj):
    """
    Write organelle stacks directly inside the Region directory (no acinus folders).
    """
    stack_dir = f"{root_dir}/stack{stackNum}"
    if not os.path.exists(stack_dir):
        os.mkdir(stack_dir)

    tifffile.imwrite(
        f"{stack_dir}/stack{stackNum}_orgs.tif",
        asinus_maxproj,
        photometric="minisblack",
    )
    print(f"Saved {stack_dir}/stack{stackNum}_orgs.tif")



# ================================================================
# Corrections – dispatcher (ONLY change here is mito vs peroxi split)
# ================================================================
def corrections(key, image, image_blur_unused):
    """
    Channel-dependent correction step.

    Parameters
    ----------
    key : str
    image : np.ndarray
        Z-stack (Z,H,W), expected uint8 for mito/peroxi path
    image_blur_unused : unused
    """
    print(f"Key received in corrections: {key}")

    if "mito" in key:
        return corrections_mito(image)

    if "peroxi" in key:
        return corrections_peroxi(image)

    # Default fallback: if routed here for any other channel, just do MIP
    return np.max(image, axis=0)




# ================================================================
# Main orchestration
# ================================================================
def main(image_path, nSlices, nStacks, channels):
    lobule_dir = image_path
    organelle_dir = image_path  # use user-provided directory directly

    # 🔴 ADD THIS BLOCK RIGHT HERE
    tifs = glob.glob(os.path.join(organelle_dir, "*.tif")) + \
           glob.glob(os.path.join(organelle_dir, "*.tiff"))

    if not tifs:
        raise FileNotFoundError(f"No TIFF files found in: {organelle_dir}")

    print(f"Found {len(tifs)} TIFF files")

    # ⬇️ existing code continues
    sample_asinus, crop_info = getSampleAsinus(organelle_dir)
    print("Sample shape:", sample_asinus.shape)

    for stackNum in range(nStacks):
        stack_dir = os.path.join(lobule_dir, f"stack{stackNum}")
        stack_images = glob.glob(f"{stack_dir}/*.tif")
        if len(stack_images) > 0:
            print(f"Stack {stackNum} already created. Skipping.")
            continue

        asinus_maxproj = getOrgStacks(
            organelle_dir, channels, nSlices, stackNum, crop_info
        )

        outputTIFs(lobule_dir, stackNum, asinus_maxproj)


# ================================================================
# Script entry point
# ================================================================
if __name__ == "__main__":
    # Single test directory (Region8)
    path = r"path_to_your_tiff"

    channels = {
        "lipid": "ch01",
        "mito": "ch00",
        "glucagon": "ch02",
        "insulin": "ch03",
        "actin": "ch04",
        "peroxi": "ch05",
    }

    n_slices = 10   # number of Z-sections per stack
    n_stacks = 5    # how many stacks (chunks) you want to process

    main(path, n_slices, n_stacks, channels)
