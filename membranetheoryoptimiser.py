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
from scipy import optimize as sc

import strake
import error
import conicalgeometrytocylindrical as geometry
import globalloads as gl
import fundamentalstressresultants as fsr
import LS3bucklingcheck as LS3

# user inputs

geoms_filename = "Sadowskietal2023-benchmarkgeometries.csv"
loads_filename = "Sadowskietal2023-benchmarkloads-LC1.csv"

# constants
f_yk = 345e6 # characteristic steel strength [N/m2]
gamma_M1 = 1.1 # EN 1993-1-6 Table 4.2
E = 210e9 # Young's Modulus [N/m2]
rho = 7850 # density [kg/m3]

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
tau_xthRd = sigma_xRd.copy()

sigma_xEd = sigma_xRd.copy()
tau_xthEd = sigma_xRd.copy()

chi_x = sigma_xRd.copy()
chi_tau = sigma_xRd.copy()

min_t = sigma_xRd.copy()

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
 
    selfWeight = np.array(fsr.p_zStresses(rho*-9.81*s.t,s.h,Z)) # selfweight as p_z

    # stresses
    N = fsr.cumulativeStresses(loads,Z,Theta,s)
    N += selfWeight

    sigma_xEd[strakeID] = np.amin(N[2,:,:])/s.t
    tau_xthEd[strakeID] = np.amin(N[1,:,:])/s.t

    # error calculations
    errors[strakeID] = error.calculateStressErrors(N,loads,s.r,tol)

    # LS3 buckling checks
    # axial
    sigma_xRd[strakeID], chi_x[strakeID] = LS3.findAxialBucklingStress(E, f_yk, fab_class, gamma_M1,s)
    checks[strakeID]["Axial check"] = LS3.checkIndividualStresses(sigma_xEd[strakeID], sigma_xRd[strakeID])

    # shear
    tau_xthRd[strakeID], chi_tau[strakeID] = LS3.findShearBucklingStress(E, f_yk, fab_class,gamma_M1,s)
    checks[strakeID]["Shear check"] = LS3.checkIndividualStresses(tau_xthEd[strakeID], tau_xthRd[strakeID])

    # interactions
    checks[strakeID]["Interaction check"] = LS3.checkStressInteractions([sigma_xEd[strakeID], 0, tau_xthEd[strakeID]],
                                                                         [sigma_xRd[strakeID],100,tau_xthRd[strakeID]],
                                                                         [chi_x[strakeID][0], 0, chi_tau[strakeID][0]])
    
    # min_t[strakeID] = sc.minimize(optimise.objectiveFunction,s.t,args=(s,loads,2,rho,E,f_yk,fab_class,gamma_M1),
                                #   bounds=sc.Bounds(1,s.r))
    
    loads["P"] += selfWeight[2,0,0]*2*math.pi*s.r # adds self weight of this can to next cans loading
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


    f.write("omega = "+str(strakes[key].omega)+"\n")
    f.write("C_x = "+str(strakes[key].C_x)+"\n")
    f.write("C_tau = "+str(strakes[key].C_tau)+"\n")
    f.write("chi_x "+str(chi_x[key][1])+"\n")
    f.write("chi_tau "+str(chi_tau[key][1])+"\n")

# f.write(str(min_t))

f.close()