#%%


os.chdir('//prfs.hhmi.org/sgrolab/mark/liver_proj/cnt_liver3/lobule_1/mito_dendra2_phall555_Lipitox_R_580_Ms_PMP70_647')

z00 = Image.open('Region 1_Merged--Z00--C00.tif')

PV1_coords = [4578,5309]
PV2_coords = [15971,2376]
PV3_coords = [18702,15464]
CV_coords = [15024,9170]

asinus1,box1 = dispCrop(z00,CV_coords,PV1_coords)
asinus2,box2 = dispCrop(z00,CV_coords,PV2_coords)
asinus3,box3 = dispCrop(z00,CV_coords,PV3_coords)



fig,ax = plt.subplots(figsize=(20,4))

ax = plt.subplot(1,4,1)
plt.imshow(z00,origin='lower')
ax.add_patch(box1)
ax.add_patch(box2)
ax.add_patch(box3)

plt.subplot(1,4,2)
plt.imshow(asinus1,origin='lower')

plt.subplot(1,4,3)
plt.imshow(asinus2,origin='lower')

plt.subplot(1,4,4)
plt.imshow(asinus3,origin='lower')