# Diagonal Crop 

import os, pickle, math, collections, numpy as np
from PIL import Image 
Image.MAX_IMAGE_PIXELS = 500000000
from scipy import ndimage
from matplotlib import patches
import tifffile
import glob

_Point = collections.namedtuple('Point', ['x', 'y'])


class Point(_Point):
    def __add__(self, p):
        return Point(self.x + p.x, self.y + p.y)

    def __sub__(self, p):
        return Point(self.x - p.x, self.y - p.y)

    def recenter(self, old_center, new_center):
        return self + (new_center - old_center)

    # http://homepages.inf.ed.ac.uk/rbf/HIPR2/rotate.htm
    def rotate(self, center, angle):
        # angle should be in radians
        x = math.cos(angle) * (self.x - center.x) - math.sin(angle) * (self.y - center.y) + center.x
        y = math.sin(angle) * (self.x - center.x) + math.cos(angle) * (self.y - center.y) + center.y
        return Point(x, y)

def getRotatedRectanglePoints(angle, base_point, height, width):
    # base_point is the upper left (ul)
    ur = Point(
        width * math.cos(angle),
        -width * math.sin(angle)
    )
    lr = Point(
        ur.x + height * math.sin(angle),
        ur.y + height * math.cos(angle)
    )
    ll = Point(
        height * math.cos(math.pi / 2 - angle),
        height * math.sin(math.pi / 2 - angle)
    )
    return tuple(base_point + pt for pt in (Point(0, 0), ur, lr, ll))

def getBounds(points):
    Bound = collections.namedtuple('Bound', ('left', 'upper', 'right', 'lower'))
    xs, ys = zip(*points)
    # left, upper, right, lower using the usual image coordinate system
    # where top-left of the image is 0, 0
    return Bound(min(xs), min(ys), max(xs), max(ys))

def roundint(values):
    return tuple(int(round(v)) for v in values)

def getBoundsCenter(bounds):
    return Point(
        (bounds.right - bounds.left) / 2 + bounds.left,
        (bounds.lower - bounds.upper) / 2 + bounds.upper
    )

def getCenter(im):
    return Point(*(d / 2 for d in im.size))

def _cropWithPoints(im, angle, points):
    bounds = getBounds(points)
    im2 = im.crop(roundint(bounds))
    bound_center = getBoundsCenter(bounds)
    crop_center = getCenter(im2)
    # in the cropped image, this is where our points are
    crop_points = [pt.recenter(bound_center, crop_center) for pt in points]
    # this is where the rotated points would end up without expansion
    rotated_points = [pt.rotate(crop_center, angle) for pt in crop_points]
    # expand is necessary so that we don't lose any part of the picture
    im3 = im2.rotate(-angle * 180 / math.pi, expand=True)
    # but, since the image has been expanded, we need to recenter
    im3_center = getCenter(im3)
    rotated_expanded_points = [pt.recenter(crop_center, im3_center) for pt in rotated_points]
    im4 = im3.crop(roundint(getBounds(rotated_expanded_points)))
    return im4

def crop(im, base, angle, height, width):
    """Return a new, cropped image.

    Args:
        im: a PIL.Image instance
        base: a (x,y) tuple for the upper left point of the cropped area
        angle: angle, in radians, for which the cropped area should be rotated
        height: height in pixels of cropped area
        width: width in pixels of cropped area
    """
    base = Point(*base)
    points = getRotatedRectanglePoints(angle, base, height, width)
    return _cropWithPoints(im, angle, points)

def getVeinCoords(lobule_dir):
    os.chdir(lobule_dir)
    with open('vein_coords.pickle','rb') as f:
        CV_coords,PV_coords = pickle.load(f)
    return np.array(CV_coords), np.array(PV_coords)

def getDAPI(
    nuclei_dir,
    nSlices,
    stackNum,
    crop_info):
    
    # get DAPI channel 
    # preallocate zstack for channel
    zstack = np.zeros((
        nSlices,
        crop_info['height'],
        crop_info['width']),
        dtype='uint8')
    
    # get files in right order
    files = glob.glob(f'{nuclei_dir}/*.tif')
    files.sort(key= lambda x: int(x.split('Z')[1].split('.')[0]))
    files.reverse()
    
    # iterate through and gather each z slice for channel
    for i in range(nSlices):
        asinus_slice = Image.open(files[nSlices*stackNum+i])
        asinus,box = dispCrop(asinus_slice, crop_info=crop_info)
        zstack[i] = asinus
    
    return ndimage.median_filter(np.max(zstack,axis=0),size=1)

def getRectangleCorner(coords,dist):
    angle = np.pi*5/4
    radius = np.sqrt(2*(.2*dist)**2)
    return [coords[0] + radius * np.cos(angle), coords[1]+radius*np.sin(angle)]

def dispCrop(img, CV_coords=None, PV_coords=None, crop_info=None):
    
    if crop_info is None:
        print('calculating crop info')
        # compute distance between portal and central vein 
        asinus_dist = np.linalg.norm(PV_coords - CV_coords)

        # compute angle between portal and central veins relative to coordinate axes 
        asinus_theta = np.arctan((CV_coords[1]-PV_coords[1])/(CV_coords[0]-PV_coords[0]))
        # if PV is to the left of CV, add pi 
        if PV_coords[0] < CV_coords[0]:
            asinus_theta += np.pi
        
        # define bottom corner of bouding rectangle 
        rectangle_coords = getRectangleCorner(CV_coords,asinus_dist)
        
        # draw box over asinus 
        box = patches.Rectangle(rectangle_coords,1.4*asinus_dist,.4*asinus_dist,angle=np.rad2deg(asinus_theta),rotation_point=(CV_coords[0],CV_coords[1]),alpha=0.5)
        
        crop_info = {
            'base': box.get_corners()[2],
            'angle': np.pi - asinus_theta,
            'height': int(0.4 * asinus_dist),
            'width': int(1.2 * asinus_dist)
        }
    # crop image according to bounding box 
    cropped_im = crop(img, crop_info['base'], crop_info['angle'], crop_info['height'], crop_info['width'])
    
    # return flipped image as array and box for image
    return np.array(cropped_im), crop_info

def getOrgStacks(
    organelle_dir,
    channels,
    nSlices,
    stackNum,
    crop_info
    ):
    
    asinus_maxproj = np.zeros((
        len(channels),
        crop_info['height'],
        crop_info['width']), 
        dtype='uint8')
    
    for i, (key, values) in enumerate(channels.items()):
        print(i)
        # get all files for channel 
        org_files = glob.glob(f'{organelle_dir}/*{values}.tif')
        org_files.sort(key= lambda x: int(x.split('Z')[1].split('-')[0]))
        org_files.reverse()

        # preallocate zstack for channel
        zstack = np.zeros((
            nSlices,
            crop_info['height'],
            crop_info['width']),
            dtype='uint8')

        # iterate through and gather each z slice for channel
        for j in range(nSlices):
            asinus_slice = Image.open(org_files[nSlices*stackNum+j])
            asinus,box = dispCrop(asinus_slice, crop_info=crop_info)
            zstack[j] = asinus
        
        asinus_maxproj[i] = ndimage.median_filter(np.max(zstack,axis=0),size=1)

    return asinus_maxproj
    
def makeAcinusDir(lobule_dir,asinusNum):
    acinus_dir = f'{lobule_dir}/acinus{asinusNum}'
    if not os.path.exists(acinus_dir):
        os.mkdir(acinus_dir)
    else:
        print('Directory already exists')
    return acinus_dir

def outputTIFs(acinus_dir,stackNum,asinusNum,DAPI_max,asinus_maxproj):
    stack_dir = f'{acinus_dir}/stack{stackNum}'
    if not os.path.exists(stack_dir):
        os.mkdir(stack_dir)
    tifffile.imwrite(f'{stack_dir}/stack{stackNum}_DAPI.tif',DAPI_max,photometric='minisblack')
    tifffile.imwrite(f'{stack_dir}/stack{stackNum}_orgs.tif',asinus_maxproj,photometric='minisblack')

def getSampleAsinus(organelle_dir,CV_coords,PV_coords,asinusNum):
    sample_image_path = glob.glob(f'{organelle_dir}/*.tif')
    sample_asinus_slice = Image.open(sample_image_path[0])
    sample_asinus, crop_info = dispCrop(
        sample_asinus_slice,
        CV_coords,
        PV_coords[asinusNum,:]
        )

    return sample_asinus, crop_info

def main(
    image_path,
    nSlices,
    nStacks,
    channels
    ):
    
    """
    Generates cropped acinus max z projections from whole lobule stiched images
    
    Args:
        image_path: path to lobule directory
        nSlices: number of slices to be assessed in max z projection 
        nStacks: number of max z projections to be cropped
    
    Returns:
        None
        Creates folder structures for each acinus and each stack under 
        image_path
    """
    
    lobule_dir = image_path
    organelle_dir = lobule_dir + '/Mito_Perox_LD_Actin'
    nuclei_dir = lobule_dir + '/DAPI'
    
    # get coordinates from file 
    CV_coords, PV_coords = getVeinCoords(lobule_dir)

    # iterate through each acinus
    for asinusNum in range(3):
        
        # create directory for acinus
        acinus_dir = makeAcinusDir(lobule_dir,asinusNum)
        
        # iterate through each z chunk
        for stackNum in range(nStacks):    

            # get sample crop for array sizing 
            sample_asinus, crop_info = getSampleAsinus(
                organelle_dir,
                CV_coords,
                PV_coords,
                asinusNum
                )

            print('ran sample')
            # get DAPI max projection 
            DAPI_max = getDAPI(
                nuclei_dir,
                nSlices,
                stackNum,
                crop_info
                )
            print('ran DAPI')

            # get all other channels 
            asinus_maxproj = getOrgStacks(
                organelle_dir,
                channels,
                nSlices,
                stackNum,
                crop_info)
            print('ran orgs')

            # output as multipage TIF files 
            outputTIFs(
                acinus_dir,
                stackNum,
                asinusNum,
                DAPI_max,
                asinus_maxproj)

if __name__ == '__main__':
    image_path = "../../../../../../../feliciano/felicianolab/For_Alex_and_Mark/Male/STV/Liv1/Lobule1"
    
    channels = {
        "actin": 'C00',
        "mito": 'C01',
        "lipid": 'C02',
        "peroxi": 'C03'
        }
    n_slices = 10
    n_stacks = 1
    main(image_path, n_slices, n_stacks, channels)