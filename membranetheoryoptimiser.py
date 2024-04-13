# membranetheoryoptimiser

"""
This script optmises the thicknesses in a tower.

@author: Cerys Morley
2024
"""
# module imports

import math
import numpy as np

import strake
import error
import conicalgeometrytocylindrical as geometry
import globalloads as gl
import fundamentalstressresultants as fsr
import LS3bucklingcheck as LS3

# user inputs

geoms_filename = "Sadowskietal2023-benchmarkgeometries.csv"
loads_filename = "Sadowskietal2023-benchmarkloads-LC1.csv"
th_n = 361
arcLength = 2*math.pi

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
loads = gl.allLoads(loads_filename,loadNames)

# geometry

strakeList, H, V = geometry.listStrakeIDs(geoms_filename)
# create strake cell in dictionary
strakes = {strakeID : strake.strake(strakeID,geoms_filename) for i, strakeID in enumerate(strakeList)}

# set up loop
theta = np.linspace(0,arcLength,th_n)
for key in strakes.keys():
    equilibriumCheck = {key : None}
    sigma_xRd = {key : None}
    tau_xthRd = {key : None}

    sigma_xEd = {key : None}
    tau_xthEd = {key : None}

    results = {key: {"Global equilibrium": None,
                     "Axial check": None,
                     "Shear check": None,
                     "Interaction check": None,
                     "Axial utilisation": None,
                     "Shear utilisation": None}}

# loop over strakes
for strakeID, s in strakes.items():
    z = np.linspace(0,s.h,100)
    Theta, Z = np.meshgrid(theta,z,indexing="ij")

    # stresses
    N = fsr.cumulativeStresses(loads,Z,Theta,s,H)
    sigma_xEd[strakeID] = np.amin(N[2,:,:])/s.t
    tau_xthEd[strakeID] = np.amin(N[1,:,:])/s.t

    # error calculations
    results[strakeID]["Global equilibrium"] = error.calculateStressErrors(N,loads,s.r,tol)

    # LS3 buckling checks
    # axial
    sigma_xRd[strakeID], chi_x = LS3.findAxialBucklingStress(E, f_yk, Q_x[fab_class], lambda_x0, chi_xh, s.h, s.r, s.t, gamma_M1)
    results[strakeID]["Axial check"], results[strakeID]["Axial utilisation"] = LS3.checkIndividualStresses(sigma_xEd[strakeID], sigma_xRd[strakeID])

    # shear
    tau_xthRd[strakeID], chi_tau = LS3.findShearBucklingStress(E, f_yk, Q_tau[fab_class], lambda_tau0, chi_tauh, beta_tau, eta_tau, s.h, s.r, s.t, gamma_M1)
    results[strakeID]["Shear check"], results[strakeID]["Shear utilisation"] = LS3.checkIndividualStresses(tau_xthEd[strakeID], tau_xthRd[strakeID])

    # interactions
    results[strakeID]["Interaction check"] = LS3.checkStressInteractions([sigma_xEd[strakeID], 0, tau_xthEd[strakeID]],
                                                                         [sigma_xRd[strakeID],100,tau_xthRd[strakeID]],
                                                                         [chi_x, 0, chi_tau])
    
# results output
for key in strakes.keys():
    print(key)
    # if all resultants are within tolerance
    if sum(results[key]["Global equilibrium"][load] for l,load in enumerate(loadNames[3:-1])) == 4:
        print("Global equilibrium satisified.")
    else:
        print("Global equilibrium violated.")
        print(results[key]["Global Equilibrium"])
    
    if results[key]["Axial check"]:
        print("Axial check passed.")
    else:
        print("Failed axial check.")
    
    print(results[key]["Axial utilisation"])

    if results[key]["Shear check"]:
        print("Shear check passed.")
    else:
        print("Failed shear check.")

    print(results[key]["Shear utilisation"])

    if results[key]["Interaction check"]:
        print("Interaction check passed.")
    else:
        print("Failed stress interaction check.")
