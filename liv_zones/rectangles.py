# Diagonal Crop 

import os, sys, pickle, numpy as np
from PIL import Image 
Image.MAX_IMAGE_PIXELS = 500000000
from matplotlib import patches
import tifffile
sys.path.append('//groups/sgro/sgrolab/mark/liver_proj/diagonal-crop')
import diagonal_crop

def getRectangleCorner(coords,dist):
    angle = np.pi*5/4
    radius = np.sqrt(2*(.2*dist)**2)
    return [coords[0] + radius * np.cos(angle), coords[1]+radius*np.sin(angle)]

def dispCrop(img,CV_coords,PV_coords):
    
    # compute distance between portal and central vein 
    asinus_dist = np.sqrt((PV_coords[0]-CV_coords[0])**2+(PV_coords[1]-CV_coords[1])**2)
    
    # compute angle between portal and central veins relative to coordinate axes 
    asinus_theta = np.arctan((CV_coords[1]-PV_coords[1])/(CV_coords[0]-PV_coords[0]))
    # if PV is to the left of CV, add pi 
    if PV_coords[0] < CV_coords[0]:
        asinus_theta += np.pi
    
    # define bottom corner of bouding rectangle 
    rectangle_coords = getRectangleCorner(CV_coords,asinus_dist)
    
    # draw box over asinus 
    box = patches.Rectangle(rectangle_coords,1.4*asinus_dist,.4*asinus_dist,angle=np.rad2deg(asinus_theta),rotation_point=(CV_coords[0],CV_coords[1]),alpha=0.5)
    
    # crop image according to bounding box 
    cropped_im = diagonal_crop.crop(img, box.get_corners()[2], np.pi-asinus_theta, .4*asinus_dist, 1.4*asinus_dist)
    
    # return flipped image as array and box for image
    return np.array(cropped_im),box


lobule_dir = '//groups/feliciano/felicianolab/For_Alex_and_Mark/Male/CNT/Liv5/Lobule3'
organelle_dir = lobule_dir + '/Mito_Perox_LD_Actin'
nuclei_dir = lobule_dir + '/DAPI'

n_zslices = 10
pixels_per_um = 22.1870

os.chdir(lobule_dir)
with open('vein_coords.pickle','rb') as f:
    CV_coords,PV_coords = pickle.load(f)

channels = ['C00','C01','C02','C03']

for asinusNum in range(3):

    # size the array
    os.chdir(organelle_dir)
    sample_asinus_slice = Image.open(os.listdir()[0])
    sample_asinus,sample_box = dispCrop(sample_asinus_slice,CV_coords,PV_coords[asinusNum])
    
    # preallocate max z projection array 
    asinus_maxproj = np.zeros([4,np.size(sample_asinus,0),np.size(sample_asinus,1)],dtype='uint8')
    
    # get DAPI channel 
    os.chdir(nuclei_dir)
    # preallocate zstack for channel
    zstack = np.zeros([n_zslices,np.size(sample_asinus,0),np.size(sample_asinus,1)],dtype='uint8')
    
    # iterate through and gather each z slice for channel
    for i in range(n_zslices):
        asinus_slice = Image.open(os.listdir()[i])
        asinus,box = dispCrop(asinus_slice,CV_coords,PV_coords[asinusNum])
        zstack[i] = asinus
    
    DAPI_max = np.max(zstack,axis=0)
    
    asinus_maxproj[1] = np.max(zstack,axis=0)
    
    # get all other channels 
    os.chdir(organelle_dir)
    for j in range(len(channels)):
        
        # get all files for channel 
        files = []
        for file in os.listdir():
            if channels[j] in file:
                files.append(file)
        
        # preallocate zstack for channel
        zstack = np.zeros([n_zslices,np.size(sample_asinus,0),np.size(sample_asinus,1)],dtype='uint8')
        
        # iterate through and gather each z slice for channel
        for i in range(n_zslices):
            asinus_slice = Image.open(files[i])
            asinus,box = dispCrop(asinus_slice,CV_coords,PV_coords[asinusNum])
            zstack[i] = asinus
        
        # index 0 for actin (c2), 1
        if j==0:
            asinus_maxproj[1] = np.max(zstack,axis=0)
        elif j==1:
            asinus_maxproj[0] = np.max(zstack,axis=0)
        elif j==2:
            asinus_maxproj[2] = np.max(zstack,axis=0)
        else:
            asinus_maxproj[3] = np.max(zstack,axis=0)
    
    # output as multipage TIF file 
    os.chdir(lobule_dir + '/z10')
    tifffile.imwrite('acinus' + str(asinusNum) + '_z' + str(n_zslices) + '_DAPI.tif',DAPI_max,photometric='minisblack')
    tifffile.imwrite('acinus' + str(asinusNum) + '_z' + str(n_zslices) + '.tif',asinus_maxproj,photometric='minisblack')








