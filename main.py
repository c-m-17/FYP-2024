"""
This script optimises the strake thicknesses in a tower.

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
import fundamentalstressresultants as fsr
import LS3bucklingcheck as LS3
import plotter
import matplotlib.pyplot as plt

# user inputs
geoms_filename : str = "CFM-geometries-01.csv"
loads_filename : str = "Sadowskietal2023-benchmarkloads-LC2.csv"

# constants
f_yk : float = 355e6 # characteristic steel strength [N/m2]
gamma_M1 : float = 1.1 # material partial factor of safety EN 1993-1-6 Table 4.2
E : float = 210e9 # Young's Modulus [N/m2]
rho : float = 7850 # density [kg/m3]

# parameters
fabClass : str = "A" # fabrication quality class (one of A, B, C)
beta_max : float = math.pi/8 # maximum apex half-angle [rad]

# loading
loadNames : list[str] = ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
loadInput : dict[str,float] = {load : geometry.importLoadMagnitude(load,loads_filename) for load in loadNames}

# geometry
strakeList, H = geometry.listStrakeIDs(geoms_filename)
strakes : dict[str | int, strake.strake] = {strakeID : strake.strake(strakeID,geoms_filename) for i, strakeID in enumerate(strakeList)}

# set up loop
theta : np.ndarray[float] = np.linspace(0,2*math.pi,361)
minimumValues : dict[str | int, sc.OptimizeResult] = {key : None for key in strakes.keys()}
loads : dict[str,float] = loadInput.copy()

# plot original geometry
fig1, (ax1, ax2) = plt.subplots(num = 1, nrows=1, ncols=2, sharey=True, layout="constrained")
plotter.plotTowerGeometry(fig1, ax1, strakes, "Initial")
ax1.set_ylabel(f"Height [m]")

def findMembraneStresses(s : strake.strake, loads : dict[str,float], rho : float) -> tuple[float,float,float]:
    """
    Calculates the critical membrane stresses throughout a strake caused by external loads and self-weight.

    s: Strake object
    loads: Value of different applied load types, with keys: [p_r, p_th, p_z, P, Q, T, M]
    rho: Material density

    returns: critical circumferential, shear and meridional stress, in that order.
    """
    theta : np.ndarray = np.linspace(0,2*math.pi,361)
    z : np.ndarray = np.linspace(0,s.h,100)
    Theta, Z = np.meshgrid(theta,z,indexing="ij")

    selfWeight : np.ndarray[float] = fsr.p_zStresses(rho*-9.81*s.t,s.h,Z)
    N : np.ndarray[float] = fsr.cumulativeStresses(loads,Z,Theta,s)
    N += selfWeight

    errors : dict[str,bool] = error.calculateStressErrors(N,loads,s.r)[0]
    for k,v in errors.items():
        if v is False:
            raise ArithmeticError("Global equilibrium is not satisfied for "+k)

    sigma_thEd : float = np.amin(N[0,:,:])/s.t
    tau_xthEd : float = np.amin(N[1,:,:])/s.t # minimum = most tensile
    sigma_xEd : float = np.amax(N[2,:,:])/s.t 

    return sigma_thEd, tau_xthEd, sigma_xEd

def objectiveFunction(x: np.ndarray, s: strake.strake, loads: dict[str,float], rho: float, E: float, f_yk: float, fabClass: str, gamma_M1=1.1) -> float:
    """
    Function which should be minimized to zero.
    Uses relative error formula to compare design stress and design resistance.

    x: Independent variable
    s: Strake object
    loads: Constant values of different applied load types named as exactly [p_r, p_th, p_z, P, Q, T, M]
    rho: Material density
    E: Young's modulus of material
    f_yk: Characteristic yield stress of material
    fabClass: Fabrication class of material, *one* of ["A","B","C"]
    gamma_M1: Material partial factor of safety, default = 1.1 (in accordance with EN 1993-1-6 Table 4.2)

    returns: Relative squared error
    """
    s.t = float(abs(x[0]))
    s.r_bot = float(abs(x[1]))
    s.update()

    sigma_thEd, tau_xthEd, sigma_xEd = findMembraneStresses(s,loads,rho)

    sigma_thRd : float = LS3.findCircumferentialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
    tau_xthRd : float = LS3.findShearBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
    sigma_xRd : float = LS3.findAxialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]

    # sign convention makes design stress match resistance sign
    relSqErr : list[float] = error.relativeSquaredError(trueValue=[sigma_thRd, tau_xthRd, sigma_xRd],
                                                        calcValue=[sigma_thEd, tau_xthEd, -sigma_xEd])
    
    f : float = min(relSqErr)

    return f

def volumeFunction(x: np.ndarray, s: strake.strake) -> float:
    """
    Returns volume of strake.

    x: Independent variable array [thickness, bottom radius]
    s: Strake object

    returns: Volume of cylindrical strake approximation.
    """

    s.t = float(abs(x[0]))
    s.r_bot = float(abs(x[1]))
    s.update()

    f : float = 2 * math.pi * ((s.r + s.t/2)**2 - (s.r - s.t/2)**2) * s.h

    return f

def resistanceCheck(x : np.ndarray, s: strake.strake, loads: dict[str,float], rho: float, E: float, f_yk: float, fabClass: str, gamma_M1=1.1) -> np.ndarray:
    """
    Calculates the relative error between the design stresses and resistances.
    To be used as a constraint for SciPy's optimize method.

    x: Independent variable
    s: Strake object
    loads: Constant values of different applied load types named as exactly [p_r, p_th, p_z, P, Q, T, M]
    rho: Material density
    E: Young's modulus of material
    f_yk: Characteristic yield stress of material
    fabClass: Fabrication class of material, *one* of ["A","B","C"]
    gamma_M1: Material partial factor of safety, default = 1.1 (in accordance with EN 1993-1-6 Table 4.2)

    returns: vector of relative errors for stress components in [th,xth,x] order,
    which must all be positive to ensure the resistance checks have passed.
    """
    s.t = float(abs(x[0]))
    s.r_bot = float(abs(x[1]))
    s.update()

    sigma_thEd, tau_xthEd, sigma_xEd = findMembraneStresses(s,loads,rho)

    sigma_thRd : float = LS3.findCircumferentialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
    tau_xthRd : float = LS3.findShearBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
    sigma_xRd : float = LS3.findAxialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]

    relErr : list[float] = error.relativeError(trueValue=[sigma_thRd, tau_xthRd, sigma_xRd],
                                               calcValue=[sigma_thEd, tau_xthEd, -sigma_xEd])

    return np.array(relErr)

def stressInteractionConstraint(x: np.ndarray, s : strake.strake, loads: dict[str,float], rho: float, E: float, f_yk: float, fabClass: str, gamma_M1=1.1) -> np.ndarray:
    """
    Checks whether the stress interaction check was passed.
    To be used as a constraint for SciPy's optimize method.

    x: Independent variable(s)
    s: Strake object
    loads: Constant values of different applied load types named as exactly [p_r, p_th, p_z, P, Q, T, M]
    rho: Material density
    E: Young's modulus of material
    f_yk: Characteristic yield stress of material
    fabClass: Fabrication class of material, *one* of ["A","B","C"]
    gamma_M1: Material partial factor of safety, default = 1.1 (in accordance with EN 1993-1-6 Table 4.2)

    returns: value of interaction check which needs to be positive to pass.
    """
    s.t = float(abs(x[0]))
    s.r_bot = float(abs(x[1]))
    s.update()

    sigma_thEd, tau_xthEd, sigma_xEd = findMembraneStresses(s,loads,rho)

    sigma_thRd, chi_theta = LS3.findCircumferentialBucklingStress(E, f_yk, fabClass, gamma_M1, s)
    tau_xthRd, chi_tau = LS3.findShearBucklingStress(E, f_yk, fabClass, gamma_M1, s)
    sigma_xRd, chi_x = LS3.findAxialBucklingStress(E, f_yk, fabClass, gamma_M1, s)

    interactionCheck : float = LS3.checkStressInteractions([sigma_thEd, tau_xthEd, sigma_xEd],
                                                          [sigma_thRd, tau_xthRd, sigma_xRd],
                                                          [chi_theta, chi_tau, chi_x])

    return np.array(interactionCheck)

jointRadius : float = list(strakes.values())[0].r_top
for strakeID, s in strakes.items():
    s.r_top = jointRadius
    s.update()

    z : np.ndarray[float] = np.linspace(0,s.h,100)
    Theta, Z = np.meshgrid(theta,z,indexing="ij")

    minimumValues[strakeID] = sc.minimize(volumeFunction,np.array([s.t,s.r_bot]),args=(s),
                                         method="SLSQP", jac="3-point", hess="3-point",
                                         constraints=[{"type":"ineq", "fun" : resistanceCheck, "args" : [s,loads,rho,E,f_yk,fabClass,gamma_M1]},
                                                      {"type":"ineq", "fun" : stressInteractionConstraint, "args" : [s,loads,rho,E,f_yk,fabClass,gamma_M1]}],
                                         bounds=[(0.002,0.1),(s.r_top, s.r_top+(s.h*math.sin(beta_max)))],
                                         options={"disp" : True})
    
    selfWeight : np.ndarray[float] = fsr.p_zStresses(rho*-9.81*s.t,s.h,Z)
    loads["P"] += selfWeight[2,0,0]*2*math.pi*s.r # adds self weight of this strake to next strake's loading
    loads["M"] += loads["Q"]*s.h # adds total moment for current strake to new baseline for next strake
    jointRadius : float = s.r_bot # saves bottom radius to be enforced as top radius of new strake

plotter.plotTowerGeometry(fig1, ax2, strakes, "Optimized")

with open(loads_filename.split(".")[0]+"-results-"+datetime.now().strftime("%y%m%d%H%M%S")+".csv","x") as f:
    f.write("Strake ID, Optimiser Success, Message, Objective Function, Minimum Thickness (mm), Minimum Bottom Radius (m), Top Radius (m),\n")
    
    for key in strakes:
        f.write(str(key)+",")
        f.write(str(minimumValues[key].get("success"))+",")
        f.write(str(minimumValues[key].get("message"))+",")
        f.write(str(minimumValues[key].get("fun"))+",")
        f.write(str(abs(minimumValues[key].get("x")[0])*1000)+",")
        f.write(str(abs(minimumValues[key].get("x")[1]))+",")
        f.write(str(strakes[key].r_top)+",\n")