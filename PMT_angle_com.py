import numpy as np
import matplotlib.pyplot as plt

"""
This script compare PMT response to beam tilt for both pencil beam and diffuser cases.
"""

### load text ###
fname = 'diff vs ang.txt'

pmt_angle, ang, diff, ang_t, diff_t = np.loadtxt(fname, skiprows=1, delimiter = ",", usecols=(0,1,2,3,4), unpack=True) #[:, :-5] 
print(diff)

### calculate Fresnel Reflection and Transmission ###
r_p = (1.33*np.cos(np.radians(pmt_angle)) - np.cos(np.arcsin(np.sin(np.radians(pmt_angle))/1.33))) / (1.33*np.cos(np.radians(pmt_angle)) + np.cos(np.arcsin(np.sin(np.radians(pmt_angle))/1.33))) #p-polarized
r_s = (np.cos(np.radians(pmt_angle)) - 1.33*np.cos(np.arcsin(np.sin(np.radians(pmt_angle))/1.33))) / (np.cos(np.radians(pmt_angle)) + 1.33*np.cos(np.arcsin(np.sin(np.radians(pmt_angle))/1.33))) #s-polarized
R_p = r_p**2
R_s = r_s**2
R_th = (R_p + R_s)/2
T_th = 1 - R_th

### calculate ratio between number of absorbed photons from beams at all angles and that at the normal ###
def ratio(a,b, angle, T):
 r1 = a/a[0]                                                                  # ratio for pencil beam case
 r2 = (b*T[0])/(b[0]*np.cos(np.radians(angle))*T)                             # ratio for diffuser case (check manual for full calculation)
 return r1, r2

r1, r2 = ratio(ang, diff, pmt_angle, T_th)
print(r1, '\n', r2)

pmt_angle = np.degrees(np.arcsin(np.sin(np.radians(pmt_angle))/1.33))
angle = np.array([pmt_angle[0], pmt_angle[1]+18.1345, pmt_angle[2]+35.3149, pmt_angle[3]+35.3149, pmt_angle[4]+18.1345, pmt_angle[5], pmt_angle[6]+18.1345, pmt_angle[7]+35.3149])           # add PMT's tilt to total beam tilting
sort = np.argsort(angle)
print("angle:", angle)
angle = np.cos(np.radians(angle))

### sorting base on total angle ###
angle = angle[sort]
diff = diff[sort]
ang = ang[sort]
ang_t = ang_t[sort]
diff_t = diff_t[sort]
r1 = r1[sort]
r2 = r2[sort]

### calculate error ###
def err(n, n_o):
   n_o_err = np.sqrt(n_o)
   n_err = np.sqrt(n)
   return (n/n_o)*np.sqrt((n_err/n)**2+(n_o_err/n_o)**2)

ang_err = err(ang, ang[0])
diff_err = err(diff, diff[0])

plt.errorbar(angle, r1, yerr=ang_err, fmt='k.', label = 'pencil beam', linestyle = '-')
plt.errorbar(angle, r2, yerr=diff_err, fmt='b.', label = 'diffuser', linestyle = '-')
plt.title("Difference between two pre-calibration methods: \n Pencil beam on a grid vs. Diffuser, injected at 50cm away from apex")
plt.xlabel("cosθ")
plt.ylabel("Ratio")
plt.legend()
plt.savefig("ang_vs_diff.png")
plt.show()