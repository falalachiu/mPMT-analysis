import numpy as np
import matplotlib.pyplot as plt

fname = 'outliers_counts.txt'

sca_x, sca_y, sca_z, sca_c, refl_x, refl_y, refl_z, refl_c = np.loadtxt(fname, skiprows = 1, delimiter = None, usecols = (0,1,2,3,4,5,6,7), unpack = True)

refl_x = np.array(refl_x)[refl_c != 0]
refl_y = np.array(refl_y)[refl_c != 0]
refl_z = np.array(refl_z)[refl_c != 0]
refl_c = np.array(refl_c)[refl_c != 0]

h = 767.57-refl_z
w = np.sqrt(refl_x**2 + refl_y**2)
angle = np.degrees(np.arctan(w/h))
rate = ((1000-refl_c)/1000)
r = 1000/rate

print(r)
print(angle)
data = np.column_stack((angle, r))
np.savetxt("angle vs reflection rate.txt", data, fmt="%f", header="Angle, ratio")

ax = plt.figure().add_subplot(111, projection='3d')
plot = ax.scatter(refl_x, refl_y, refl_z, c=r, s=3)
ax.set_xlabel('x position (mm)')
ax.set_ylabel('y position (mm)')
ax.set_zlabel('z position (mm)')
ax.set_title('Photon reflecting ratio vs. beam injection angle')
ax.view_init(elev=90, azim=-90)
plt.colorbar(plot, label='Number of reflected photons', shrink=0.6)
plt.savefig('Reflected_angle.png', dpi=300)
plt.show()
plt.close()

plt.errorbar(angle, r,  fmt = 'k.')
plt.title("Photon reflecting ratio vs. beam injection angle")
plt.xlabel("Angle")
plt.ylabel("Adjusted rate")
plt.savefig('Reflected_angle_rate.png')
plt.show()