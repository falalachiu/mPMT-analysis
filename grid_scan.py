## run.mac editor 
## comment out unused part based on sections
import numpy as np
import sympy as sp

######### Simple angular check (towards origin) ##########
# it is set to do movement in y direction now with shift so it will end up at the origin all time
# angle = 20                                                                   
# y = np.sin(np.radians(angle))*(300/np.cos(np.radians(angle)))      # y-shift

# with open("build/run.mac", "w") as f:
#     f.write(f"/mygenerator/SetX 0 \n")                             # set x position
#     f.write(f"/mygenerator/SetY {y} \n")                           # set y position
#     f.write(f"/mygenerator/SetZ 767.57 \n")                        # set z position
#     f.write(f"/mygenerator/SetPZenithAngle {angle} \n")            # set zenith angle (angle from -z axis)
#     f.write(f"/mygenerator/SetPAzimuthAngle -90 \n")               # set azimuth angle (0 is towards +x direction)
#     f.write("/run/beamOn 1000000\n\n")                             # set how many beams of photons is injected in this setup

########## pencil beam scan ##########
# angle data (written in excel but put here for convenience in editing)
# (0 -> +1) 11.765081851751518, (0 -> +2) 19.76624897040672
# (+1 -> 0) 12.0896077610188, (+1 -> +2) 9.11847778484164, (+1 -> -1)22.654529969829795, (+1 -> -2) 29.545728336654577

# zen = 0                                                                            # injection angle (manual or use above data to edit)
# zen2 = np.arcsin(np.sin(np.radians(zen))/1.33)                                     # beam tilt in water when there's separate water layer (just Snell's Law calculation)
# def shift(x):                                                                      # calculation of x-shift
#     ### when have separate water layer
#     # ex_depth = 200 + 347 - np.sqrt(347**2 - x**2)                                # extra depth in z-direction (compared to apex of mPMT)
#     # xpos = x - 300 * np.tan(np.radians(zen)) - ex_depth * np.tan(zen2)     

#     ### when fully in water
#     xpos = x - 500 * np.tan(np.radians(zen)) 
#     return xpos
# xshift_min = shift(-45)                                                            # shifted minimum x value (data in excel file)
# xshift_max = shift(45)+1                                                           # shifted maximum x value (data in excel file)
# x_val = np.arange(xshift_min, xshift_max, 0.1)                                     # set x range 
# y_val = np.arange(-45, 46, 0.1)                                                    # set y range

# with open("build/run.mac", "w") as f:
#     for x in x_val:
#         for y in y_val:
#                 f.write(f"/mygenerator/SetX {x}\n")
#                 f.write(f"/mygenerator/SetY {y}\n")
#                 f.write(f"/mygenerator/SetZ 767.57 \n")
#                 # f.write(f"/mygenerator/SetDiscRad 1 \n")                         # set beam radius
#                 # f.write(f"/mygenerator/SetSpread 5 \n")                          # set beam angle spread (*could be a replacement for diffusive scan but this will record data of different directions in 1 list, need new logic to read data)
#                 # f.write(f"/mygenerator/SetAngle 0 \n")                           # not setup properly (don't use it)
#                 f.write(f"/mygenerator/SetPZenithAngle {zen} \n")
#                 f.write(f"/mygenerator/SetPAzimuthAngle 0 \n")
#                 f.write("/run/beamOn 10\n\n")

# ### check-ins ###
# print(np.degrees(zen2), xshift_min, xshift_max)
# # print(np.tan(np.radians(zen))*400, np.tan(zen2)*100)
# print(len(x_val), len(y_val))

######### Fibonacci point in 3D ##########

# def fibonacci_pt(n_pt, max_zen, R):
#     #set golden angle and lin z space for uniform distribution
#     phi_golden = np.pi * (3 - np.sqrt(5))
#     z_min = np.cos(np.radians(max_zen))

#     #calculate points
#     pts = []
#     for i in range(n_pt):
#         z_norm = 1 - (i / float(n_pt - 1)) * (1 - z_min)
#         zen_rad = np.arccos(z_norm)
#         azi_rad = phi_golden * i

#         x = R * np.sin(zen_rad) * np.cos(azi_rad)
#         y = R * np.sin(zen_rad) * np.sin(azi_rad)
#         z = R * np.cos(zen_rad) - 191.43

#         pts.append({
#             'x': x,
#             'y': y,
#             'z': z,
#             'theta': np.degrees(zen_rad),
#             'psi': (np.degrees(azi_rad) % 360)
#         })
#     return pts

# R = 691.43
# num_pt = 200000
# max_zenith = 73

# scan_pts = fibonacci_pt(num_pt, max_zenith, R)

# with open("build/run.mac", "w") as f:
#     for pt in scan_pts: 
#         zenith = pt['theta']
#         azimuth = pt['psi'] + 180
#         if azimuth > 180:
#             azimuth -= 360
#         elif azimuth < -180:
#             azimuth += 360

#         if pt['z'] <= 0: 
#             continue

#         f.write(f"/mygenerator/SetX {pt['x']}\n")
#         f.write(f"/mygenerator/SetY {pt['y']}\n")
#         f.write(f"/mygenerator/SetZ {pt['z']}\n")
#         f.write(f"/mygenerator/SetPZenithAngle {zenith} \n")
#         f.write(f"/mygenerator/SetPAzimuthAngle {azimuth} \n")
#         f.write("/run/beamOn 1000\n\n")

#     print(f"Total scan points: {len(scan_pts)}")

######### Diffusive scan ##########
def shift(x):                                                                                             # Calculation for x range for target PMT
    a = sp.Symbol('a')
    ex_depth = 200 + 347 - np.sqrt(347**2 - (x-100)**2)                                                   # extra depth; use (x-n) for diffuser at different position 
    eqn = (1.33**2 - 1) * (a * (x - a))**2 + 1.33**2 * (300 * (x - a))**2 - (ex_depth * a)**2
    soln = sp.solve(eqn, a)
    for sol in soln:
        if sol.is_real and sol < x:
            return float(sol)
    return sol

min_x = shift(1)                                                              # shifted x-range (check excel for data), when including vertical beam, do not put 0 for min_x, use 1 and edit other parts below to prevent calculation error
max_x = shift(350)
print(min_x, max_x)
min_height = 500+347-np.sqrt(347**2-min_x**2)
max_height = 500+347-np.sqrt(347**2-max_x**2)
zen_min = np.degrees(np.arctan(min_x/300))                                        # calculation zenith angle range, replace min_x with 0 when include vertical beam 
zen_max = np.degrees(np.arctan(max_x/300))

scan_pt = round((500000/(255**2 * np.pi))*(350**2*np.pi))                     # set number of beams to inject (usually manual input,current calculation is for full )
print(scan_pt)

psi = np.degrees(np.arctan(50/175))                                           # set azimuth range (change denominator for calculation, denominator = horizontal distance between diffuser and target PMT centre)
min_psi = -psi                                                                # when including vertical beam, it must be a 360º scan, use min_psi = -180, max_psi = 180
max_psi = psi

## checking ###
zen = (zen_min + zen_max)/2
print(zen)
print(min_height, max_height, zen_min, zen_max, psi) 

with open("build/run.mac", "w") as f:
    for i in range(scan_pt):
        zen = np.random.uniform(zen_min, zen_max)
        psi = np.random.uniform(min_psi, max_psi)

        f.write(f"/mygenerator/SetX 100\n")                                                   # set x position (diffuser positon)
        f.write(f"/mygenerator/SetY 0\n")
        f.write(f"/mygenerator/SetZ 767.57 \n")
        f.write(f"/mygenerator/SetPZenithAngle {zen} \n")              
        f.write(f"/mygenerator/SetPAzimuthAngle {psi} \n")
        f.write("/run/beamOn 100\n\n")