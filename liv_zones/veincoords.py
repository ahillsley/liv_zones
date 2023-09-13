# creating pickle files for vein coordinates

import os, pickle

os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv1/Lobule1')
CV_coords = [7630,10532]
PV_coords =[[5537,2410],[13948,3961],[3017,18517]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)
    
os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv1/Lobule2')
CV_coords = [9523,6929]
PV_coords =[[2511,2176],[7193,15114],[19092,10861]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)
    
os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv2/Lobule1')
CV_coords = [10820,10644]
PV_coords =[[15061,2268],[19991,6501],[2939,14323]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)

os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv2/Lobule2')
CV_coords = [10463,3403]
PV_coords =[[2273,2864],[5770,9985],[17981,5667]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)

os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv2/Lobule3')
CV_coords = [9043,2683]
PV_coords =[[2613,9073],[14240,9547],[18447,2791]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)
    
os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv3/Lobule1')
CV_coords = [14906,9117]
PV_coords =[[4513,5404],[15840,2413],[18744,15415]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)

os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv3/Lobule2')
CV_coords = [13500,8930]
PV_coords =[[3814,4331],[4237,12304],[17468,16241]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)

os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv3/Lobule3')
CV_coords = [10679,12176]
PV_coords =[[2055,11588],[13951,3022],[19665,15857]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)
    
os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv5/Lobule1')
CV_coords = [8956,10495]
PV_coords =[[2169,4228],[13550,2587],[6644,17987]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)
    
os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv5/Lobule2')
CV_coords = [8564,13759]
PV_coords =[[5112,3239],[20552,6879],[2759,23497]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)

os.chdir('//prfs.hhmi.org/felicianolab/For_Alex_and_Mark/Male/CNT/Liv5/Lobule3')
CV_coords = [13124,5200]
PV_coords =[[2070,3413],[10042,13481],[22398,6646]]
with open('vein_coords.pickle','wb') as f:
    pickle.dump([CV_coords,PV_coords],f,pickle.HIGHEST_PROTOCOL)


