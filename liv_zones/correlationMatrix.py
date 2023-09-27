# Mito Analysis

import os
import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt 
from matplotlib import image as img 
import scipy 

os.chdir('//prfs.hhmi.org/sgrolab/mark/liver_proj/sampledata')

mitoData = pd.read_csv('mitochondria_properties.csv')
ldData = pd.read_csv('lipid_dropplet_properties.csv')
cellData = pd.read_csv('average_properties_per_cell.csv')
labelcells = img.imread('labled_cells.png')

cellDataClean = cellData.drop(38)
properties = cellDataClean.columns[[1,5,6,7,8,9,10,
                                   11,13,15,16,17,18,19,
                                   21,22,24,25,26,27,28,
                                   30,31,32,33,34,35,39,
                                   40,41,48,49,
                                   54]]


rvals = np.zeros([len(properties),len(properties)])

for i in range(len(properties)):
    for j in range(0,i+1):
        fitData = scipy.stats.linregress(cellDataClean[properties[i]],cellDataClean[properties[j]])
        
        if fitData[3] < 0.01:
            rvals[i,j] = fitData[2]
        else:
            rvals[i,j] = 'NaN'   


plt.figure(figsize=(12,10))
plt.title('$R^2$ correlation (p<0.01)',fontsize=24)
cmap = plt.imshow(rvals,vmin=-1,vmax=1,cmap='coolwarm')
plt.colorbar(cmap)
plt.xticks(np.linspace(0,len(properties)-1,len(properties)),properties,rotation=45,ha='right')
plt.yticks(np.linspace(0,len(properties)-1,len(properties)),properties)
plt.tight_layout()
    



