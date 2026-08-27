import numpy as np
import matplotlib.pyplot as plt
"""
This script compares the absorption probabilities of runs with different wavelength setups. 
"""
### Set 1 (less dense wavelengths range with scan in different PMTs) ###
# fname = 'wavelength_response.txt'

# wv, pmt0, pmt1, pmt2 = np.loadtxt(fname, skiprows=1, delimiter = ",", usecols=(0,1,2,3), unpack=True) 

# pmt0_err = np.sqrt(pmt0)
# pmt1_err = np.sqrt(pmt1)
# pmt2_err = np.sqrt(pmt2)

# angle = np.array([0, 11.765081851751518+18.1345, 19.76624897040672+35.3149])                          # add PMT tilt angle to get the full beam tilt compared to beam perpendicular to mPMT apex
# print(angle)

# plt.errorbar(wv, pmt0, yerr=pmt0_err, fmt='k.', label = 'PMT0 (0º)', linestyle = '-')
# plt.errorbar(wv, pmt1, yerr=pmt1_err, fmt='b.', label = 'PMT1 (29.89958185º)', linestyle = '-')
# plt.errorbar(wv, pmt2, yerr=pmt2_err, fmt='g.', label = 'PMT2 (55.08114897º)', linestyle = '-')

### Set 2 (finer wavelengths scan all with same PMT) ###
fname = 'wv_range_response.txt'

wv, pmt = np.loadtxt(fname, skiprows=1, delimiter = ",", usecols=(0,1), unpack=True)

pmt_err = np.sqrt(pmt)
plt.errorbar(wv, pmt, yerr=pmt_err, fmt='k.', label = 'PMT - vertical beam', linestyle = '-')

#######################################################################################################################

plt.title("Relationship between PMT responses and photon energy")
plt.xlabel("Wavelength (nm)")
plt.ylabel("absorption probability")
plt.legend()
plt.savefig("wv_res.png")
plt.show()