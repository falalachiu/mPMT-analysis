import numpy as np
import matplotlib.pyplot as plt

# file 1: Reflection on PMTs; unpolarized
# file 2: Reflection on spots; with local polarizations

fname = 'Reflection_results.txt'
angle, ph_total, ph_reflected, ph_t_diff, ph_rfl_diff = np.loadtxt(fname, skiprows=1, delimiter = ",", usecols=(0,1,2,3,4), unpack=True)
R_test_ang = ph_reflected/ph_total
R_test_diff = ph_rfl_diff/ph_t_diff
error = np.sqrt(ph_reflected)/ph_total
err = np.sqrt(ph_rfl_diff)/ph_t_diff

# fname2 = 'Reflection_local.txt'
# angle, s_ph, p_ph, un_ph = np.loadtxt(fname2, skiprows=1, delimiter = ",", usecols=(0,1,2,3), unpack=True)
# cos = np.cos(np.radians(angle))
# R_test_s = s_ph/1000000
# R_test_p = p_ph/1000000
# R_test = un_ph/1000000
# R_test_com = (s_ph/1000000 + p_ph/1000000)/2
# print("R_test_s:", R_test_s, "\n R_test_p:", R_test_p, "\n R_test_com:", R_test_com)
# serr = np.sqrt(s_ph)/1000000
# perr = np.sqrt(p_ph)/1000000
# unerr = np.sqrt(un_ph)/1000000
# comerr = np.sqrt(s_ph+p_ph)/1000000/2

### combined theoretical calculation ###

r_p = (1.34*np.cos(np.radians(angle)) - np.cos(np.arcsin(np.sin(np.radians(angle))/1.34))) / (1.34*np.cos(np.radians(angle)) + np.cos(np.arcsin(np.sin(np.radians(angle))/1.34))) #p-polarized
r_s = (np.cos(np.radians(angle)) - 1.34*np.cos(np.arcsin(np.sin(np.radians(angle))/1.34))) / (np.cos(np.radians(angle)) + 1.34*np.cos(np.arcsin(np.sin(np.radians(angle))/1.34))) #s-polarized
R_p = r_p**2
R_s = r_s**2
print("R_p:", R_p, "\n R_s:", R_s)
# R_theta = (np.sin(np.radians(angle)))**2*R_p + (np.cos(np.radians(angle)))**2*R_s
R_th = (R_p + R_s)/2

### plotting ###

#1
plt.errorbar(angle, R_test_ang, yerr=error, fmt = 'k.', label = 'simulated (pencil beam; unpolarized)')
plt.errorbar(angle, R_test_diff, yerr=err, fmt = 'b.', label = 'simulated (diffuser; unpolarized)')
plt.scatter(angle, R_th, label = 'theoretical unpolarized', marker = 'x')

#2
# plt.errorbar(angle, R_test_s, yerr=serr, fmt = 'k.', label = 'simulated (local s-direction)')
# plt.errorbar(angle, R_test_p, yerr=perr, fmt = 'b.', label = 'simulated (local p-direction)')
# plt.errorbar(angle, R_test, yerr=unerr, fmt = 'r.', label = 'simulated (local unpolarized)')
# # plt.errorbar(angle, R_test_com, yerr=comerr, fmt = 'g.', label = 'simulated (local combined)')
# plt.scatter(angle, R_s, label = 'theoretical s-polarized', marker = 'x')
# plt.scatter(angle, R_p, label = 'theoretical p-polarized', marker = 'x')
# plt.scatter(angle, R_th, label = 'theoretical unpolarized', marker = 'x')
# # plt.scatter(angle, R_theta, label = 'theoretical angularly-polarized', marker = 'x')

plt.title("Theoretical vs. Simulated reflection ratio of photons at water surface \n for Fresnel reflection at different angles")
plt.xlabel("angle θ (º)")
plt.ylabel("Reflection ratio")
plt.legend()
plt.savefig("Reflection_ratio.png")
