from xml.parsers.expat import errors
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Ellipse
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from scipy import ndimage
import uproot
import numpy as np
import matplotlib.pyplot as plt
import os  
import itertools

""" 
This script compares data from 2 runs from raw data, no data extraction from extract_2D.py is needed. 
"""

### calls out 2 data files ###
fname = os.path.join(
    '/home', 'fchiu', 'projects', 'def-blairt2k', 'fchiu', 'outputs', 'output_mPMT492.root'
)

fname2 = os.path.join(
    '/home', 'fchiu', 'projects', 'def-blairt2k', 'fchiu', 'outputs', 'output_mPMT494.root' 
)

### extract absorption data from both files (check inline comments in extract_2D.py for details) ### 
def extract_absorption_ratios(filename, filename2):
    def grab_data(data):
        ratios = []
        errors = []
        ph_total = []
        ph_absorbed = []

        pos_x = []
        pos_y = []
        pos_z = []
    
        for i in data.keys():
            if i.startswith("Photons_") and "Master" not in i:                                                      
                all_x = data[i]["PosX"].array()                             
                all_y = data[i]["PosY"].array()
                all_z = data[i]["PosZ"].array()
                step = data[i]["Step_Number"].array()
                id = data[i]["Photon_ID"].array()

                int_mask = (step == 2)
                land_x = all_x[int_mask]                                                                                                
                land_y = all_y[int_mask]
                land_z = all_z[int_mask]
                counts_x = Counter(land_x)
                counts_y = Counter(land_y)
                ref_x = counts_x.most_common(1)[0][0]
                ref_y = counts_y.most_common(1)[0][0]
                valid_mask = (np.abs(land_x - ref_x) < 1) & (np.abs(land_y - ref_y) < 1) & (150 < land_z) & (land_z < 267.57)
                sca_mask = (np.abs(land_x - ref_x) >= 1) & (np.abs(land_y - ref_y) >= 1) & (land_z > 155.57) & (land_z < 467.57)
                valid_x = land_x[valid_mask]
                valid_y = land_y[valid_mask]
                valid_z = land_z[valid_mask]

                pos_x.append(np.mean(valid_x))
                pos_y.append(np.mean(valid_y))
                pos_z.append(np.mean(valid_z))                

                n_photons = len(np.unique(data[i]["Photon_ID"].array())) # total number of photons (by reading the number marking of photons in root)
                n_absorbed = np.sum(data[i]["Step_Status"]. array() == 2) # number of absorbed photons
                ph_total.append(n_photons)
                ph_absorbed.append(n_absorbed)

                p = (n_absorbed / n_photons)    # absorption ratio
                ratios.append(p) # add absorption ratio to the ratio data list

                err = np.sqrt(p)/n_photons
                errors.append(err) #add error to the error data list
        return np.array(ratios), np.array(errors), np.array(ph_total), np.array(ph_absorbed), np.array(pos_x), np.array(pos_y), np.array(pos_z)


    data1 = uproot.open(filename)
    data2 = uproot.open(filename2)

    ratio, error, ph_total, ph_absorbed, pos_x, pos_y, pos_z = grab_data(data1)
    ratio2, error2, ph_total2, ph_absorbed2, pos_x2, pos_y2, pos_z2 = grab_data(data2)

    ### calculate difference between 2 data ###
    ratio_diff = ratio2/(ratio+0.00001)                                 # target data over control data (add 0.00001 to denominator to prevent error at 0 reading points)
    # ratio_diff = ratio - ratio2                                       # compare by looking at difference in count

    return ratio, error, ph_total, ph_absorbed, pos_x, pos_y, pos_z, ratio2, error2, ph_total2, ph_absorbed2, pos_x2, pos_y2, pos_z2, ratio_diff

### save absorption data ###
ratio, error, ph_total, ph_absorbed, pos_x, pos_y, pos_z, ratio2, error2, ph_total2, ph_absorbed2, pos_x2, pos_y2, pos_z2, ratio_diff = extract_absorption_ratios(fname, fname2)    

save_txt = np.column_stack((pos_x, pos_y, ratio_diff, ratio, ratio2))
np.savetxt("compared.txt", save_txt, fmt="%f", header="x position, y position, difference in absorption ratio, ratio 1, ratio 2")

### plot compared results ### 
ax = plt.figure().add_subplot()
ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(pos_x, pos_y, pos_z, c=ratio_diff, cmap='viridis', vmin=0, vmax=1, s=5, alpha=1.0, edgecolors='none')
ax.set_aspect('equal')
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)') 
ax.set_zlabel('z position (mm)')
ax.set_title('Difference in Absorption Ratio (360/340)')
ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)
ax.set_zlim(150, 270)

ax.set_box_aspect([1, 1, 1.2])
plt.colorbar(plot, label='Difference in absorption Ratio relative to 340nm plot', shrink=0.6)
ax.view_init(elev=90, azim=-90)
plt.savefig('abs_diff.png', dpi=300)
plt.show()
plt.close()