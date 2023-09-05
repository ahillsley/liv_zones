# Diagonal Crop 

import os, numpy as np
from PIL import Image 
Image.MAX_IMAGE_PIXELS = 500000000
from matplotlib import pyplot as plt 
from matplotlib import patches

os.chdir('//prfs.hhmi.org/sgrolab/mark/liver_proj/diagonal-crop')
import diagonal_crop

def getRectangleCorner(coords,dist):
    angle = np.pi*5/4
    radius = np.sqrt(2*(.2*dist)**2)
    return [coords[0] + radius * np.cos(angle), coords[1]+radius*np.sin(angle)]

def dispCrop(CV_coords,PV_coords):
    asinus_dist = np.sqrt((PV_coords[0]-CV_coords[0])**2+(PV_coords[1]-CV_coords[1])**2)
    asinus_theta = np.arctan((CV_coords[1]-PV_coords[1])/(CV_coords[0]-PV_coords[0]))+np.pi
    rectangle_coords = getRectangleCorner(CV_coords,asinus_dist)
    
    box = patches.Rectangle(rectangle_coords,1.4*asinus_dist,.4*asinus_dist,angle=np.rad2deg(asinus_theta),rotation_point=(CV_coords[0],CV_coords[1]),alpha=0.5)

    return diagonal_crop.crop(z00, box.get_corners()[2], np.pi-asinus_theta, .4*asinus_dist, 1.4*asinus_dist)
    


os.chdir('//prfs.hhmi.org/sgrolab/mark/liver_proj/cnt_liver3/lobule_1/mito_dendra2_phall555_Lipitox_R_580_Ms_PMP70_647')

pixels_per_um = 22.1870


z00 = Image.open('Region 1_Merged--Z00--C00.tif')

PV1_coords = [4578,5309]
PV2_coords = [15971,2376]
CV_coords = [15024,9170]

asinus1 = dispCrop(CV_coords,PV1_coords)
asinus2 = dispCrop(CV_coords,PV2_coords)

plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.imshow(asinus1)

plt.subplot(2,1,2)
plt.imshow(asinus2)


#%%

plt.figure(figsize=(12,8))








#%%

fig,ax = plt.subplots(figsize=(8,6))


plt.imshow(z00,origin='lower')
plt.scatter(PV1_coords[0],PV1_coords[1],color='white',s=10)
plt.scatter(CV_coords[0],CV_coords[1],color='white',s=10)
plt.scatter(corner_coords[0],corner_coords[1],color='yellow',s=10)



ax.add_patch(box)
# ax.add_patch(patches.Rectangle(corner_coords,1.4*A1_dist,.4*A1_dist,angle=-angle/np.pi*180,alpha=0.5))

# plt.scatter(corner_coords2[0],corner_coords2[1],color='r',s=10)




