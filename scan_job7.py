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

fname = os.path.join(
    '/home', 'fchiu', 'projects', 'def-blairt2k', 'fchiu', 'outputs', 'output_mPMT393.root'
)
def extract_absorption_ratios(filename):
    data = uproot.open(filename)

    int_x = []
    int_y = []
    int_z = []

    ratios = []
    errors = []
    ph_total = []
    ph_absorbed = []
    #ph_reflect = []
    pos_x = []
    pos_y = []
    pos_z = []

    end_x = []
    end_y = []
    end_z = []  

    out_x = []
    out_y = []
    out_z = []

    refl_x, refl_y, refl_z = [], [], []
    sca_x, sca_y, sca_z = [], [], []
    w_x, w_y, w_z = [], [], []

    for i in data.keys():
    #     if i.startswith("Photons_Master"):
    #         int_x = data[i]["PosX_Initial"].array()
    #         int_y = data[i]["PosY_Initial"].array()
    #         int_z = data[i]["PosZ_Initial"].array()

        if i.startswith("Photons_") and "Master" not in i:
            all_x = data[i]["PosX"].array()
            all_y = data[i]["PosY"].array()
            all_z = data[i]["PosZ"].array()
            step = data[i]["Step_Number"].array()
            id = data[i]["Photon_ID"].array()

            surf_x = np.mean(all_x[step == 1])
            surf_y = np.mean(all_y[step == 1])
            surf_z = np.mean(all_z[step == 1])

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
            valid_z = land_z[valid_mask]  # use x coor to to match datapoint
            ex_x = land_x[~valid_mask]
            ex_y = land_y[~valid_mask]
            ex_z = land_z[~valid_mask]   # use x coor to to match datapoint
            scatter_x, scatter_y, scatter_z = land_x[sca_mask], land_y[sca_mask], land_z[sca_mask]
            

            stat = data[i]["Step_Status"].array()
            refl_mask = (stat == 1) & (step == 1)
            reflected_x, reflected_y, reflected_z = all_x[refl_mask], all_y[refl_mask], all_z[refl_mask]

            id1 = id[step == 1]
            w_id2 = id[step == 2][land_z == 467.57]
            x1, y1, z1 = all_x[step == 1], all_y[step == 1], all_z[step == 1]
            weird_mask = np.isin(id1, w_id2) & (z1 == 467.57)
            weird_x, weird_y, weird_z = x1[weird_mask], y1[weird_mask], z1[weird_mask]

            int_x.append(surf_x)
            int_y.append(surf_y)
            int_z.append(surf_z)
            pos_x.append(np.mean(valid_x))
            pos_y.append(np.mean(valid_y))
            pos_z.append(np.mean(valid_z))
            out_x.append(ex_x)
            out_y.append(ex_y)
            out_z.append(ex_z)
            refl_x.append(reflected_x)
            refl_y.append(reflected_y)
            refl_z.append(reflected_z)
            sca_x.append(scatter_x)
            sca_y.append(scatter_y)
            sca_z.append(scatter_z)
            w_x.append(weird_x)
            w_y.append(weird_y)
            w_z.append(weird_z)

            # print("landing position before filtering:", land_x, land_y, land_z)
            print("reference landing position:", ref_x)
            # print("valid landing position:", valid_x, valid_y, valid_z)
            # print("extra landing position:", ex_x, ex_y, ex_z)

            fin_x = data[i]["PosX"].array()[-1]
            fin_y = data[i]["PosY"].array()[-1]
            fin_z = data[i]["PosZ"].array()[-1]
            if data[i]["Step_Status"].array()[-1] == 2: # if the last step is absorption, add the position to the end position list
                end_x.append(fin_x)
                end_y.append(fin_y)
                end_z.append(fin_z) 

            n_photons = len(np.unique(data[i]["Photon_ID"].array())) # total number of photons (by reading the number marking of photons in root)
            n_absorbed = np.sum(data[i]["Step_Status"]. array() == 2) # number of absorbed photons

            ph_total.append(n_photons)
            ph_absorbed.append(n_absorbed)
            #ph_reflect.append(n_reflect)

            p = n_absorbed / n_photons if n_absorbed > 0 else 0   # absorption ratio
            ratios.append(p) # add absorption ratio to the ratio data list

            err = np.sqrt(p * (1 - p) / n_photons) if n_absorbed > 0 else 0 #error
            errors.append(err) #add error to the error data list

    end_x = np.array(end_x)
    end_y = np.array(end_y)
    end_z = np.array(end_z) 
    # print(end_x, end_y, end_z)
    pos_x = np.array(pos_x)
    pos_y = np.array(pos_y)
    pos_z = np.array(pos_z)

    dis_x = np.abs(pos_x - int_x)
    dis_y = np.abs(pos_y - int_y)
    dis_z = np.abs(pos_z - int_z)
    dis = np.sqrt(dis_x**2 + dis_y**2)
    cos = np.cos(np.arctan(dis/dis_z))

    ratios = np.array(ratios)
    ratios = ratios*cos
    ph_absorbed= np.array(ph_absorbed)*cos
    ph_total = np.array(ph_total)*cos
    # ph_absorbed = np.array(ph_absorbed)
    # ph_total = np.array(ph_total)
    errors = np.array(errors)
    
    # ph_reflect = np.array(ph_reflect)
    # print(ph_reflect)

    return ratios, errors, ph_total, ph_absorbed, pos_x, pos_y, pos_z, out_x, out_y, out_z, cos, end_x, end_y, end_z, refl_x, refl_y, refl_z, sca_x, sca_y, sca_z, w_x, w_y, w_z#, ph_reflect

#save absorption data 
# ratios, errors, ph_total, ph_absorbed = extract_absorption_ratios(fname)    
ratios, errors, ph_total, ph_absorbed, pos_x, pos_y, pos_z, out_x, out_y, out_z, cos, end_x, end_y, end_z, refl_x, refl_y, refl_z, sca_x, sca_y, sca_z, w_x, w_y, w_z = extract_absorption_ratios(fname)   

####### 2D Cartisian ########
# pos_x = np.arange(500, 601, 1)
# pos_y = np.arange(-45, 46, 1)
# ratios = np.reshape(ratios, (len(pos_x), len(pos_y)))
# errors = np.reshape(errors, (len(pos_x), len(pos_y)))
# ph_absorbed = np.reshape(ph_absorbed, (len(pos_x), len(pos_y)))
# ph_total = np.reshape(ph_total, (len(pos_x), len(pos_y)))

# print(len(pos_x), len(pos_y))

# x_centers = 0.5*(pos_x[1:] + pos_x[:-1])
# y_centers = 0.5*(pos_y[1:] + pos_y[:-1])

# print(fname)

# for x in x_centers:
#     for y in y_centers:
#         # simulation... 
#         pass

# X_coor, Y_coor = np.meshgrid(pos_x, pos_y, indexing='ij')

# xmin = 75
# xmax = 130
# corner = ((X_coor < xmin) & (Y_coor < -40)) | ((X_coor < xmin) & (Y_coor > 40)) | ((X_coor > xmax) & (Y_coor < -40)) | ((X_coor > xmax) & (Y_coor > 40))
# ph_absorbed[corner]=0
# ratios[corner]=0
# errors[corner]=0

####### 2D Spherical ########

# pos_x = []
# pos_y = []
# zenith = np.arange(0, 74, 1)
# azimuth = np.arange(-180, 181, 3)
# Z, A = np.meshgrid(np.radians(zenith), np.radians(azimuth), indexing='ij')
# R = 691.43
# pos_x = R * np.sin(Z)*np.cos(A)
# pos_y = R * np.sin(Z)*np.sin(A)
# full_r = np.zeros(pos_x.shape)
# full_err = np.zeros(pos_x.shape)
# abs = np.zeros(pos_x.shape)
# tot = np.zeros(pos_x.shape)

# counter = 0
# for i, theta in enumerate(zenith):
#     for j, psi in enumerate(azimuth): 
        
#         x = R*np.sin(np.radians(theta))*np.cos(np.radians(psi))
#         y = R*np.sin(np.radians(theta))*np.sin(np.radians(psi))

#         if theta == 0:
#             full_r[i,j] = ratios[0]
#             full_err[i,j] = errors[0]
#             abs[i,j] = ph_absorbed[0]
#             tot[i,j] = ph_total[0]

#         else:
#             if counter < len(ratios):
#                 full_r[i,j] = ratios[counter]
#                 full_err[i,j] = errors[counter]
#                 abs[i,j] = ph_absorbed[counter]
#                 tot[i,j] = ph_total[counter]
#                 counter += 1

# ratios = np.reshape(ratios, (len(zenith), len(azimuth)))
# errors = np.reshape(errors, (len(zenith), len(azimuth)))
# ph_absorbed = np.reshape(ph_absorbed, (len(zenith), len(azimuth)))
# ph_total = np.reshape(ph_total, (len(zenith), len(azimuth)))
# ph_reflect = np.reshape(ph_reflect, (len(zenith), len(azimuth)))


# save_data = np.column_stack((np.degrees(Z).flatten(), np.degrees(A).flatten(), pos_x.flatten(), pos_y.flatten(), tot.flatten(), abs.flatten(), full_r.flatten(),full_err.flatten()))
# np.savetxt("absorption_data.txt", save_data, fmt="%f", header="zenith angle (º), azimuth angle (º), x position (mm), y position (mm), total number of photons, number of absorbed photons, absorption ratios , errors")
# save_data = np.column_stack((X_coor.flatten(), Y_coor.flatten(), ph_total.flatten(), ph_absorbed.flatten(), ratios.flatten(),errors.flatten()))
# np.savetxt("absorption_data.txt", save_data, fmt="%f", header="x position (mm), y position (mm), total number of photons, number of absorbed photons, absorption ratios , errors")

# plt.figure()
# #plt.colorbar(plt.tripcolor(pltx, plty, pltr, shading = 'flat', vmax=1, vmin=0, cmap='viridis')) #
# plt.colorbar(plt.pcolormesh(pos_x, pos_y, ratios.T, shading= 'nearest', vmax=1, vmin=0, cmap='viridis'))
# plt.xlabel('x position (mm)')
# plt.ylabel('y position (mm)')
# plt.savefig('absorption_ratios_2D_new.png')
# #plt.gca().set_aspect('equal') #
# plt.show() 

# x = pos_x.flatten()
# y = pos_y.flatten()
# r = full_r.flatten()

# pt = np.vstack((x,y)).T
# _,plt_pt = np.unique(np.round(pt, 8), axis = 0, return_index = True)
# pltx, plty, pltr = x[plt_pt], y[plt_pt], r[plt_pt]
# grid_x, grid_y = np.mgrid[pltx.min():pltx.max():1000j,plty.min():plty.max():1000j]

# grid_r = griddata((pltx, plty), pltr, (grid_x, grid_y), method = 'nearest')
# disc = np.sqrt(grid_x**2+ grid_y**2)
# max_rad = np.max(np.sqrt(pltx**2+plty**2))
# grid_r[disc > max_rad] = np.nan

#######  3D plot #######

header="x position (mm), y position (mm), z position (mm), cosine of the angle, total number of photons, number of absorbed photons, absorption ratios , errors, outliers: (x,y,z)"
with open("absorption_data.txt", "w") as f:
    f.write(header + "\n")
    for i in range(len(pos_x)):
        sim_row = (f"{pos_x[i]}, {pos_y[i]}, {pos_z[i]}, {cos[i]}, {ph_total[i]}, {ph_absorbed[i]}, {ratios[i]}, {errors[i]}")
        com_sets = [(out_x, out_y, out_z), (refl_x, refl_y, refl_z), (sca_x, sca_y, sca_z)]
        for spx, spy, spz in com_sets:
            com_row = []
            if i < len(spx) and spx[i] is not None and hasattr(spx[i], "__iter__"):
                for x,y,z in zip(spx[i], spy[i], spz[i]):
                    com_row.append(f"({x}, {y}, {z})")
            if com_row:
                sim_row += ", " + "; ".join(com_row)
        f.write(sim_row + "\n")

header = "x position (mm), y position (mm), z position (mm), reflected positions(x,y,z), scattered positions(x,y,z)"
with open("outliers.txt", "w") as f:
    f.write(header + "\n")
    for i in range(len(pos_x)):
        pos_row = (f"{pos_x[i]}, {pos_y[i]}, {pos_z[i]}")
        out_sets = [(refl_x, refl_y, refl_z), (sca_x, sca_y, sca_z)]
        for spx, spy, spz in out_sets:
            out_row = []
            if i < len(spx) and spx[i] is not None and hasattr(spx[i], "__iter__"):
                for x,y,z in zip(spx[i], spy[i], spz[i]):
                    out_row.append(f"({x}, {y}, {z})")
            if out_row:
                pos_row += ", " + "; ".join(out_row)
        f.write(pos_row + "\n")            

pos_x = np.array(pos_x)
pos_y = np.array(pos_y)
pos_z = np.array(pos_z)
ratios = np.array(ratios)
print(pos_x, pos_y, pos_z)

end_x = np.array(end_x)
end_y = np.array(end_y)
end_z = np.array(end_z)

# R_lim = 691.43
# disc = np.sqrt(pos_x**2 + pos_y**2 + (pos_z+191.43)**2)
# mask = np.abs(disc - R_lim) < 20
# mask = pos_z > 155.57
# pos_x = pos_x[mask]
# pos_y = pos_y[mask]
# pos_z = pos_z[mask]
# ratios = ratios[mask]
# ph_absorbed = ph_absorbed[mask]
# ph_total = ph_total[mask]
# errors = errors[mask]

# xmin = 150
# xmax = 190
# corner = ((pos_x < xmin) & (pos_y < -40)) | ((pos_x < xmin) & (pos_y > 40)) | ((pos_x > xmax) & (pos_y < -40)) | ((pos_x > xmax) & (pos_y > 40))
# ph_absorbed[corner]=0
# ph_total[corner]=0
# ratios[corner]=0
# errors[corner]=0

# valid_area = Ellipse(xy = (0, 0), width = 90, height = 90, edgecolor = 'red', fill = False, alpha = 0.6, lw = 2)

ax = plt.figure().add_subplot()
ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(pos_x, pos_y, pos_z, c=ratios, cmap='viridis', vmin=0, vmax=1, s=5, alpha=1.0, edgecolors='none')
# plot = ax.scatter(pos_x, pos_y, c=ratios, cmap='viridis', vmin=0, vmax=1, s=5, alpha=1.0, edgecolors='none')
# ax.add_patch(valid_area)
ax.set_aspect('equal')
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)') 
ax.set_zlabel('z position (mm)')
ax.set_title('Absorption Ratios')

ax.set_xlim(-200, -140)
ax.set_ylim(-50, 50)
ax.set_zlim(150, 270)

ax.set_box_aspect([0.6, 1, 1.2])
plt.colorbar(plot, label='Absorption Ratio', shrink=0.6)
ax.view_init(elev=90, azim=-90)
plt.savefig('absorption_ratios_3D_top-down_7.png', dpi=300)
plt.show()
plt.close()

ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(pos_x, pos_y, pos_z, c=ratios, cmap='viridis', vmin=0, vmax=1, s=5, alpha=1.0, edgecolors='none')
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)')
ax.set_zlabel('z position (mm)')
ax.set_title('Absorption Ratios in 3D')
ax.set_xlim(-200, -140)
ax.set_ylim(-50, 50)
ax.set_zlim(150, 270)

ax.set_box_aspect([0.6, 1, 1.2])
plt.colorbar(plot, label='Absorption Ratio', shrink=0.6)
# ax.view_init(elev=90, azim=-90)
plt.savefig('absorption_ratios_3D_7.png', dpi=300)
plt.show()

tilt = 0
print("Total number of absorbed photons:", np.sum(ph_absorbed)*np.cos(np.radians(tilt)))
print("Integrated absorption ratio:", (np.sum(ratios)/len(ratios))*np.cos(np.radians(tilt)))

# def valid_counter(ph_absorbed):
#     mask = (ph_absorbed > 10) & (np.sqrt(((pos_x-100)/45)**2 + (pos_y/45)**2) < 1)
#     return np.sum(mask)
 
mask = np.sqrt((((pos_x+167.5)/25)**2 + (pos_y/42.5)**2)) < 1
valid_pt = np.sum(mask)
print("Total pixels with absorption:", valid_pt)
print("Total number of photons absorbed in valid pixels:", np.sum(ph_absorbed[mask])*np.cos(np.radians(tilt)))
print("Integrated absorption ratio:", (np.sum(ratios[mask])/valid_pt)*np.cos(np.radians(tilt)))

end_x = np.array(end_x)
end_y = np.array(end_y)
end_z = np.array(end_z)

ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(end_x, end_y, end_z, s=3)
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)')
ax.set_zlabel('z position (mm)')
ax.set_title('Absorption spots top view')
ax.view_init(elev=90, azim=-90)
plt.savefig('absorption_spots_3D_7.png', dpi=300)
plt.show()

#####################################################################################

out_x_mask = [out_x[i] for i in range(len(out_x)) if mask[i]]
out_y_mask = [out_y[i] for i in range(len(out_y)) if mask[i]]
out_z_mask = [out_z[i] for i in range(len(out_z)) if mask[i]]   

refl_x_mask = [refl_x[i] for i in range(len(refl_x)) if mask[i]]
refl_y_mask = [refl_y[i] for i in range(len(refl_y)) if mask[i]]
refl_z_mask = [refl_z[i] for i in range(len(refl_z)) if mask[i]]   

sca_x_mask = [sca_x[i] for i in range(len(sca_x)) if mask[i]]
sca_y_mask = [sca_y[i] for i in range(len(sca_y)) if mask[i]]
sca_z_mask = [sca_z[i] for i in range(len(sca_z)) if mask[i]]

w_x_mask = [w_x[i] for i in range(len(w_x)) if mask[i]]
w_y_mask = [w_y[i] for i in range(len(w_y)) if mask[i]]
w_z_mask = [w_z[i] for i in range(len(w_z)) if mask[i]]

out_x = np.array(list(itertools.chain.from_iterable(out_x)))
out_y = np.array(list(itertools.chain.from_iterable(out_y)))
out_z = np.array(list(itertools.chain.from_iterable(out_z)))

refl_x = np.array(list(itertools.chain.from_iterable(refl_x_mask)))
refl_y = np.array(list(itertools.chain.from_iterable(refl_y_mask)))
refl_z = np.array(list(itertools.chain.from_iterable(refl_z_mask)))

sca_x = np.array(list(itertools.chain.from_iterable(sca_x_mask)))
sca_y = np.array(list(itertools.chain.from_iterable(sca_y_mask)))
sca_z = np.array(list(itertools.chain.from_iterable(sca_z_mask)))

w_x = np.array(list(itertools.chain.from_iterable(w_x_mask)))
w_y = np.array(list(itertools.chain.from_iterable(w_y_mask)))
w_z = np.array(list(itertools.chain.from_iterable(w_z_mask)))

print("Total number of outliers:", len(out_x))
print("Total number of reflected photons:", len(refl_x))
print("Total number of scattered photons:", len(sca_x))
print("Total number of weird photons:", len(w_x_mask))

sca = np.column_stack((sca_x, sca_y, sca_z))
refl = np.column_stack((refl_x, refl_y, refl_z))
w = np.column_stack((w_x, w_y, w_z))

def dupli_count(set):
    values, counts = np.unique(set, axis = 0, return_counts = True)
    return values[:,0], values[:,1], values[:,2], counts

sca_x_v, sca_y_v, sca_z_v, sca_c = dupli_count(sca)
refl_x_v, refl_y_v, refl_z_v, refl_c = dupli_count(refl)
w_x_v, w_y_v, w_z_v, w_c = dupli_count(w)
# print("Reflected photon x direction values:", refl_x_v)
# print("Reflected photon x direction counts:", refl_c)

ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(sca_x_v, sca_y_v, sca_z_v, c=sca_c, s=3)
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)')
ax.set_zlabel('z position (mm)')
ax.set_title('Scattered photons spots in 3D')
ax.view_init(elev=90, azim=-90)
plt.colorbar(plot, label='Number of scattered photons', shrink=0.6)
plt.savefig('scattered_spots_3D_7.png', dpi=300)
plt.show()
plt.close()

ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(refl_x_v, refl_y_v, refl_z_v, c=refl_c, s=3)
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)')
ax.set_zlabel('z position (mm)')
ax.set_title('Reflected photon spots top view')
ax.view_init(elev=90, azim=-90)
plt.colorbar(plot, label='Number of reflected photons', shrink=0.6)
plt.savefig('Reflected_spots_top-down_7.png', dpi=300)
plt.show()

ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(w_x_v, w_y_v, w_z_v, c=w_c, s=3)
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)')
ax.set_zlabel('z position (mm)')
ax.set_title('Weird photon spots top view')
ax.view_init(elev=90, azim=-90)
plt.colorbar(plot, label='Number of reflected photons', shrink=0.6)
plt.savefig('Weird_spots_top-down_7.png', dpi=300)
plt.show()

if len(sca_x) < len(refl_x):
    sca_x = np.pad(sca_x, (0,len(refl_x)-len(sca_x)), mode = 'constant', constant_values=0)
    sca_y = np.pad(sca_y, (0,len(refl_y)-len(sca_y)), mode = 'constant', constant_values=0)
    sca_z = np.pad(sca_z, (0,len(refl_z)-len(sca_z)), mode = 'constant', constant_values=0)

save_data = np.column_stack((refl_x, refl_y, refl_z, sca_x, sca_y, sca_z))
np.savetxt("outliers_list_7.txt", save_data, fmt="%f", header="reflection photons position (x,y,z), scattered photon position (x,y,z)")

if len(sca_x_v) < len(refl_x_v):
    sca_x_v = np.pad(sca_x_v, (0,len(refl_x_v)-len(sca_x_v)), mode = 'constant', constant_values=0)
    sca_y_v = np.pad(sca_y_v, (0,len(refl_y_v)-len(sca_y_v)), mode = 'constant', constant_values=0)
    sca_z_v = np.pad(sca_z_v, (0,len(refl_z_v)-len(sca_z_v)), mode = 'constant', constant_values=0)
    sca_c = np.pad(sca_c, (0,len(refl_c)-len(sca_c)), mode = 'constant', constant_values=0)
 
counter = np.column_stack((sca_x_v, sca_y_v, sca_z_v, sca_c, refl_x_v, refl_y_v, refl_z_v, refl_c))
np.savetxt("outliers_counts_7.txt", counter, fmt="%f", header="scattered photons position (x,y,z, counts), reflected photon position (x,y,z, counts)")

angle = np.degrees(np.arctan(np.sqrt(w_x_v**2+w_y_v**2)/(300**2)))
weird_data = np.column_stack((angle, w_x_v, w_y_v, w_z_v, w_c))
np.savetxt("weird photons_7.txt", weird_data, fmt="%f", header="weird photons injection angle tilt, landing positions (x,y,z) , count")