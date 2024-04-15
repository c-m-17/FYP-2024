# membranetheoryoptimiser

"""
This script optmises the thicknesses in a tower.

@author: Cerys Morley
2024
"""
# module imports

import math
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

import strake
import error
import conicalgeometrytocylindrical as geometry
import globalloads as gl
import fundamentalstressresultants as fsr
import LS3bucklingcheck as LS3

# user inputs

geoms_filename = "Sadowskietal2023-benchmarkgeometries.csv"
loads_filename = "Sadowskietal2023-benchmarkloads-LC2.csv"

# constants
f_yk = 345 # characteristic steel strength [N/mm2]
gamma_M1 = 1.1 # EN 1993-1-6 Table 4.2
E = 210e3 # Young's Modulus [N/mm2]
rho = 7850e-9 # density [kg/mm3]
g = -9.81 # gravity [m/s2]

# axial parameters
Q_x = {"Class A":40, "Class B":25, "Class C": 16} # EN 1993-1-6 Table D.1
chi_xh = 1.10 # EN 1993-1-6 D.19
lambda_x0 = 0.10 # EN 1993-1-6 D.10

# shear parameters
Q_tau = {"Class A":40, "Class B":25, "Class C": 16} # EN 1993-1-6 Table D.7
chi_tauh = 1.0 # EN 1993-1-6 D.53
lambda_tau0 = 0.40 # EN 1993-1-6 D.50
beta_tau = 0.60 # EN 1993-1-6 D.51
eta_tau = 1.0 # EN 1993-1-6 D.52

# parameters
tol = 1e-3
fab_class = "Class A"

# loading

loadNames = ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
loadinput = gl.allLoads(loads_filename,loadNames)

# geometry

strakeList, H, V = geometry.listStrakeIDs(geoms_filename)
# create strake cell in dictionary
strakes = {strakeID : strake.strake(strakeID,geoms_filename) for i, strakeID in enumerate(strakeList)}

# set up loop
theta = np.linspace(0,2*math.pi,361)
sigma_xRd = {key : None for key in strakes.keys()}
tau_xthRd = {key : None for key in strakes.keys()}

sigma_xEd = {key : None for key in strakes.keys()}
tau_xthEd = {key : None for key in strakes.keys()}

# {strakeID : {checktype : list[boolean, values]}}
checks = {key: {"Axial check": [None]*2,
                "Shear check": [None]*2,
                "Interaction check": [None]*2} for key in strakes.keys()}

# {strakeID : [{loadtype : boolean}, {loadtype : error}]}
errors = {key: [None]*2 for key in strakes.keys()}

loads = loadinput.copy()
# loop over strakes
for strakeID, s in strakes.items():
    z = np.linspace(0,s.h,100)
    Theta, Z = np.meshgrid(theta,z,indexing="ij")
 
    selfWeight = np.array(fsr.p_zStresses(rho*g*s.t,s.h,Z)) # selfweight as p_z

    # stresses
    N = fsr.cumulativeStresses(loads,Z,Theta,s)
    N += selfWeight

    sigma_xEd[strakeID] = np.amin(N[2,:,:])/s.t
    tau_xthEd[strakeID] = np.amin(N[1,:,:])/s.t

    # error calculations
    errors[strakeID] = error.calculateStressErrors(N,loads,s.r,tol)

    # LS3 buckling checks
    # axial
    sigma_xRd[strakeID], chi_x = LS3.findAxialBucklingStress(E, f_yk, Q_x[fab_class], lambda_x0, chi_xh, s.h, s.r, s.t, gamma_M1)
    checks[strakeID]["Axial check"] = LS3.checkIndividualStresses(sigma_xEd[strakeID], sigma_xRd[strakeID])

    # shear
    tau_xthRd[strakeID], chi_tau = LS3.findShearBucklingStress(E, f_yk, Q_tau[fab_class], lambda_tau0, chi_tauh, beta_tau, eta_tau, s.h, s.r, s.t, gamma_M1)
    checks[strakeID]["Shear check"] = LS3.checkIndividualStresses(tau_xthEd[strakeID], tau_xthRd[strakeID])

    # interactions
    checks[strakeID]["Interaction check"] = LS3.checkStressInteractions([sigma_xEd[strakeID], 0, tau_xthEd[strakeID]],
                                                                         [sigma_xRd[strakeID],100,tau_xthRd[strakeID]],
                                                                         [chi_x, 0, chi_tau])
    
    loads["P"] += np.sum(selfWeight) # adds self weight of can to next cans loading
    loads["M"] += loads["Q"]*s.h # adds moment for current strake to new baseline for next strake

# checks output

f = open(loads_filename.split(".")[0]+"-checks-"+datetime.now().strftime("%y%m%d%H%M%S")+".txt","x")
for key in strakes.keys():
    f.write("\n"+str(key)+"\n")
    # if all resultants are within tolerance
    if sum(errors[key][0].values()) == len(errors[key][0].keys()):
        f.write("Global equilibrium satisified.\n")
    else:
        f.write("Global equilibrium violated.\n")
    
    f.write(str(errors[key][1])+"\n") # print error values
    
    for check in checks[key].keys():
        error.printBoolCheck(f,check,checks[key][check][0],checks[key][check][1])

f.close()