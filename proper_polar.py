import numpy as np
import matplotlib.pyplot as plt

"""
This script looks at the change detection efficiency as injection angle changes for different polarization of light. 
Although named s and p polarization, but it is just the naming, it does not refer to vectors perpendicular to the plane, but rather 2 vectors orthogonal to beam. 
"""
### load text file ###
fname = 'proper_polar.txt'
pmt_angle, spol, ppol = np.loadtxt(fname, skiprows=1, delimiter = ",", usecols=(0,1,2), unpack=True)

unpol = (spol+ppol)/2                    # calculate unpolarized case using the average of the 2 orthogonal polarizations

### recalculate beam angle when there's a separate water layer ### 
# pmt_angle = np.degrees(np.arcsin(np.sin(np.radians(pmt_angle))/1.33))
# angle = np.array([pmt_angle[0], pmt_angle[1]+18.1345, pmt_angle[2]+35.3149, pmt_angle[3]+35.3149, pmt_angle[4]+18.1345, pmt_angle[5], pmt_angle[6]+18.1345, pmt_angle[7]+35.3149])

### simply use the injection angles when using full water volume and at same PMT (no angle change) ###
angle = pmt_angle
print("angle:", angle)
angle = np.cos(np.radians(angle))

### Calculate polarization errors ###
s_err = np.sqrt(spol)
p_err = np.sqrt(ppol)
un_err = np.sqrt((np.sqrt(spol)**2 + np.sqrt(ppol)**2))

plt.errorbar(angle, spol, yerr=s_err, fmt='k.', label = 's-polarized')
plt.errorbar(angle, ppol, yerr=p_err, fmt='b.', label = 'p-polarized')
plt.errorbar(angle, unpol, yerr=un_err, fmt='g.', label = 'unpolarized')
# plt.plot(angle, ry, '-', label = 'y-polarized', marker = '.')
# plt.plot(angle, rx, '-', label = 'x-polarized', marker = '.')
# plt.plot(angle, r30, '-', label = '30º down-polarized', marker = '.')
# plt.plot(angle, r60, '-', label = '60º down-polarized', marker = '.')
plt.title("Diffrence between PMT response on s-polarized light and p-polarized \n light injected towards different PMTs in pencil beams")
plt.xlabel("cosθ")
plt.ylabel("Number of absorbed photon")
plt.legend()
plt.savefig("s_vs_p(pencil).png")
plt.show()