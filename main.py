"""
This script optimises the volume of strakes in a wind turbine support tower.

Cerys Morley
2024
"""
# package imports
import math
import numpy as np
import time
from scipy import optimize as sc
import pandas as pd
import argparse

# module imports
from strake import strake
import LS3bucklingcheck as LS3

# functions
def listStrakeIDs(filename: str) -> pd.Series:
    """
    Extracts strake IDs from file. Calculates tower height and volume.

    filename: File name (in "inputs" folder)

    returns: strake IDs, tower height[m]
    """
    df : pd.DataFrame = pd.read_csv(f"inputs\{filename}")
    StrakeIDlist : pd.Series = df["ID"]

    return StrakeIDlist

def importLoadMagnitude(loadType: str, filename: str) -> float:
    """
    Extracts the magnitude of a specific load type from a file.
    File must have "Load" & "Magnitude" columns.

    loadType: One of [p_r, p_th, p_z, P, T, Q, M]
    filename: File name (in "inputs" folder)

    returns: Global constant value of that load.
    """
    
    loads : pd.DataFrame = pd.read_csv(f"inputs\{filename}")
    A : pd.Series = loads["Magnitude"].loc[loads["Load"] == loadType]
    magnitude : float = A.tolist()[0]

    return magnitude

def p_zStresses(p_z: float, h: float, Z: np.ndarray[float,float]) -> np.ndarray:
    """
    Calculates contribution to Ns from applied tangential pressure.

    p_z: Applied tangential pressure
    h: Strake height
    Z: Local z-coordinate values corresponding to strake, in interval [0,strake height]

    returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
    correspond to [component, theta-coord, z-coord].
    """
    dN_th = np.zeros(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = p_z*(h-Z) # in z-coords
    return np.array([dN_th, dN_zth, dN_z])

def cumulativeStresses(loads: dict[str,float], Z: np.ndarray[float,float], Theta: np.ndarray[float,float], s: strake) -> np.ndarray:
    """
    Sums the contributions to membrane stress resultants from all load types.

    loads: constant values of different applied load types named as exactly [p_r, p_th, p_z, P, Q, T, M]
    Z: Local z-coordinate values corresponding to strake, in interval [0,strake height]
    Theta: Local th-coordinate values corresponding to strake, in interval [0,2*pi]
    s: strake object
    
    returns: cumulative N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
    correspond to [N component, theta-coord, z-coord].
    """
    # nested functions
    def p_rStresses(p_r: float, R: float, Z: np.ndarray[float,float]) -> np.ndarray:
        """
        Calculates contribution to Ns from applied normal pressure.

        p_r: Applied normal pressure (internal pressure = positive)
        R: Strake radius
        Z: Local z-coordinate values corresponding to strake, in interval [0,strake height]

        returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
        correspond to [component, theta-coord, z-coord].
        """
        dN_th = p_r*R*np.ones(np.shape(Z))
        dN_zth = np.zeros(np.shape(Z))
        dN_z = np.zeros(np.shape(Z))
        return np.array([dN_th, dN_zth, dN_z])

    def p_thStresses(p_th: float, Z: np.ndarray[float,float]) -> np.ndarray:
        """
        Calculates contribution to Ns from applied circumferential pressure.

        p_th: Applied circumferential pressure
        Z: Local z-coordinate values corresponding to strake, in interval [0,strake height]

        returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
        correspond to [component, theta-coord, z-coord].
        """
        dN_th = np.zeros(np.shape(Z))
        dN_zth = -p_th*Z
        dN_z = np.zeros(np.shape(Z))
        return np.array([dN_th, dN_zth, dN_z])

    def PStresses(P: float, R: float, Z: np.ndarray[float,float]) -> np.ndarray:
        """
        Calculates contribution to Ns from applied vertical force.

        P: Vertical point load, applied at centre-top of tower.
        R: Strake radius
        Z: Local z-coordinate values corresponding to strake, in interval [0,strake height]

        returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
        correspond to [component, theta-coord, z-coord].
        """
        dN_th = np.zeros(np.shape(Z))
        dN_zth = np.zeros(np.shape(Z))
        dN_z = P*np.ones(np.shape(Z))/(2*math.pi*R)
        return np.array([dN_th, dN_zth, dN_z])

    def TStresses(T: float, R: float, Theta: np.ndarray[float,float]) -> np.ndarray:
        """
        Calculates contribution to Ns from applied torque.

        T: Applied global torque (around z-axis).
        R: Strake radius
        Theta: Local th-coordinate values corresponding to strake, in interval [0,2*pi]

        returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
        correspond to [component, theta-coord, z-coord].
        """
        dN_th = np.zeros(np.shape(Theta))
        dN_zth = T*np.ones(np.shape(Theta))/(2*math.pi*R**2)
        dN_z = np.zeros(np.shape(Theta))
        return np.array([dN_th, dN_zth, dN_z])

    def MStresses(M: float, R: float, Theta: np.ndarray[float,float]) -> np.ndarray:
        """
        Calculates contribution to Ns from applied global bending moment.

        M: Applied global moment, around radial-axis, perpendicular to shear force Q direction.
        R: Strake radius
        Theta: Local th-coordinate values corresponding to strake, in interval [0,2*pi]

        returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
        correspond to [component, theta-coord, z-coord].
        """
        dN_th = np.zeros(np.shape(Theta))
        dN_zth = np.zeros(np.shape(Theta)) # *applied* moment gradient = 0. gradient from Q included there.
        dN_z = M*np.cos(Theta)/(math.pi*R**2) # cos(th)
        return np.array([dN_th, dN_zth, dN_z])

    def QStresses(Q: float, R: float, leverArm: np.ndarray[float,float], Theta: np.ndarray[float,float]) -> np.ndarray:
        """
        Calculates contribution to Ns from applied shear force.

        Q: Applied shear force, at centre-top of tower, in radial direction
        R: Strake radius
        leverArm: Distance between top of tower and z-coordinate values
        Theta: Local th-coordinate values corresponding to strake, in interval [0,2*pi]

        returns: N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
        correspond to [component, theta-coord, z-coord].
        """
        dN_th = np.zeros(np.shape(Theta))
        dN_zth = Q*np.sin(Theta)/(math.pi*R) # sin(th)
        dN_z = -Q*np.cos(Theta)*leverArm/(math.pi*R*R) # cos(th)
        return np.array([dN_th, dN_zth, dN_z])

    # initialise array
    N : np.ndarray[float] = np.zeros((3,*np.shape(Z)))

    # add stresses
    N += p_rStresses(loads["p_r"],s.r,Z)
    N += p_thStresses(loads["p_th"],Z)
    N += p_zStresses(loads["p_z"],s.h,Z)
    N += PStresses(loads["P"],s.r,Z)
    N += QStresses(loads["Q"],s.r,s.h-Z,Theta)
    N += TStresses(loads["T"],s.r,Theta)
    N += MStresses(loads["M"],s.r,Theta)

    return N

def findMembraneStresses(s : strake, loads : dict[str,float], rho : float) -> tuple[float,float,float]:
    """
    Calculates the critical membrane stresses throughout a given strake.

    s: Strake object
    loads: Value of different applied external load types, including all keys: ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
    rho: Material density (constant)

    returns: Critical circumferential, shear and meridional stress, in that order.
    """
    # set up coordinate arrays
    theta : np.ndarray = np.linspace(0,2*math.pi,361)
    z : np.ndarray = np.linspace(0,s.h,100)
    Theta, Z = np.meshgrid(theta,z,indexing="ij")

    # calculate & superimpose membrane stresses from loading
    N : np.ndarray[float] = cumulativeStresses(loads,Z,Theta,s)
    selfWeight : np.ndarray[float] = p_zStresses(rho*-9.81*s.t,s.h,Z)
    N += selfWeight

    # check global equilibrium
    errors : dict[str,bool] = calculateStressErrors(N,loads,s.r)[0]
    for k,v in errors.items():
        if v is False:
            raise ArithmeticError("Global equilibrium is not satisfied for "+k)

    # find critical value for each component (most tensile)
    sigma_thEd : float = np.amin(N[0,:,:]).tolist()/s.t
    tau_xthEd : float = np.amin(N[1,:,:]).tolist()/s.t  
    sigma_xEd : float = np.amax(N[2,:,:]).tolist()/s.t 

    return sigma_thEd, tau_xthEd, sigma_xEd

def volumeFunction(x: np.ndarray, s: strake) -> float:
    """
    Returns volume of strake after updating strake with new geometry.

    x: Independent variable array [thickness, bottom radius]
    s: Strake object

    returns: Volume of cylindrical strake approximation.
    """
    # update strake object with new geometry
    s.t = float(x[0])
    s.update(r_bot=float(x[1]))

    # calculate volume
    V : float = 2 * math.pi * ((s.r + s.t/2)**2 - (s.r - s.t/2)**2) * s.h

    return V

def resistanceCheck(x : np.ndarray, s: strake, loads: dict[str,float], rho: float, E: float, f_yk: float, fabClass: str, gamma_M1=1.1) -> np.ndarray:
    """
    Calculates the relative error between the design stresses and resistances.
    To be used as a constraint for SciPy's optimize method.

    x: Independent variables: [thickness, bottom radius]
    s: Strake object
    loads: Values of different applied external load types, including all keys: ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
    rho: Material density
    E: Young's modulus of material
    f_yk: Characteristic yield stress of material
    fabClass: Fabrication class of material, *one* of ["A","B","C"]
    gamma_M1: Material partial factor of safety, default = 1.1 (in accordance with EN 1993-1-6 Table 4.2)

    returns: vector of relative errors for stress components in [xth,x] order,
    which must all be positive to ensure the resistance checks have passed.
    """
    # update strake object with new geometry
    s.t = float(x[0])
    s.update(r_bot=float(x[1]))

    # calculate design stresses
    sigma_thEd, tau_xthEd, sigma_xEd = findMembraneStresses(s,loads,rho)

    # calculate design resistances
    sigma_thRd : float = LS3.findCircumferentialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
    tau_xthRd : float = LS3.findShearBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
    sigma_xRd : float = LS3.findAxialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]

    # calculate relative errors
    relErr : list[float] = relativeError(trueValue=[sigma_thRd, tau_xthRd, sigma_xRd],
                                         calcValue=[sigma_thEd, tau_xthEd, -sigma_xEd])

    return np.array(relErr)

def stressInteractionConstraint(x: np.ndarray, s : strake, loads: dict[str,float], rho: float, E: float, f_yk: float, fabClass: str, gamma_M1=1.1) -> np.ndarray:
    """
    Checks whether the stress interaction check was passed.
    To be used as a constraint for SciPy's optimize method.

    x: Independent variables [thickness, bottom radius]
    s: Strake object
    loads: Values of different applied external load types, including all keys: ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
    rho: Material density
    E: Young's modulus of material
    f_yk: Characteristic yield stress of material
    fabClass: Fabrication class of material, *one* of ["A","B","C"]
    gamma_M1: Material partial factor of safety, default = 1.1 (in accordance with EN 1993-1-6 Table 4.2)

    returns: value of interaction check which must be positive to pass.
    """
    # update strake object with new geometry
    s.t = float(x[0])
    s.update(r_bot=float(x[1]))

    # calculate design stresses
    sigma_thEd, tau_xthEd, sigma_xEd = findMembraneStresses(s,loads,rho)

    # calculate design resistances
    sigma_thRd, chi_theta = LS3.findCircumferentialBucklingStress(E, f_yk, fabClass, gamma_M1, s)
    tau_xthRd, chi_tau = LS3.findShearBucklingStress(E, f_yk, fabClass, gamma_M1, s)
    sigma_xRd, chi_x = LS3.findAxialBucklingStress(E, f_yk, fabClass, gamma_M1, s)

    # calculate value of buckling interaction check
    interactionCheck : float = LS3.checkStressInteractions([sigma_thEd, tau_xthEd, sigma_xEd],
                                                          [sigma_thRd, tau_xthRd, sigma_xRd],
                                                          [chi_theta, chi_tau, chi_x])

    return np.array(interactionCheck)

def relativeError(*,trueValue : float|int|list, calcValue : float|int|list) -> float | list[float]:
    """
    Calculates the relative error between two values.

    trueValue: exact value, which is also used to normalise error (unless it's equal to zero).
    calcValue: calculated value being checked

    returns: normalised decimal error value
    """
    # remove int as dtype going forward
    if type(trueValue) is int:
        trueValue = float(trueValue)
    
    if type(calcValue) is int:
        calcValue = float(calcValue)

    # checks list/float type mismatch
    if type(trueValue) is not type(calcValue):
        raise TypeError("Input arguments must have same dtype.")
    
    try:
        err : float = (trueValue - calcValue)/abs(trueValue)
    except TypeError:
        if len(trueValue) != len(calcValue):
            raise IndexError("Lists must be of the same length.")
        
        err : list = []
        for i in range(len(trueValue)):
            try:
                err.append((trueValue[i] - calcValue[i])/abs(trueValue[i]))
            except ZeroDivisionError:
                err.append((trueValue[i] - calcValue[i]))
    except ZeroDivisionError:
        err : float = (trueValue - calcValue)

    return err

def calculateStressErrors(N : np.ndarray[float], loads : dict[str,float], r : float, tol : float=1e-10) -> tuple[dict[str,bool],dict[str,float]]:
    """
    Checks global equilibrium by checking theoretical membrane stress resultants against
    calculated ones.

    N: calculated membrane stress resultants in [theta, shear, z] components in that order.
    loads: value of different applied load types in order: [p_r, p_th, p_z, P, Q, T, M]
    r: strake radius
    tol: error tolerance (default = 1e-10)
    
    returns: Global equilibrium passed? True/False, Relative error value
    """
    def relativeSquaredError(*,trueValue : float | int | list, calcValue : float | int | list) -> float | list:
        """
        Calculates the relative squared error between two values.

        trueValue: exact value, which is also used to normalise error (unless it's equal to zero).
        calcValue: calculated value being checked

        returns: normalised decimal error value
        """
        # remove int as dtype going forward
        if type(trueValue) is int:
            trueValue = float(trueValue)
        
        if type(calcValue) is int:
            calcValue = float(calcValue)

        # checks list/float type mismatch
        if type(trueValue) is not type(calcValue):
            raise TypeError("Input arguments must have same dtype.")
        
        try:
            err : float = ((trueValue - calcValue)**2)/(trueValue**2)
        except TypeError:
            if len(trueValue) != len(calcValue):
                raise IndexError("Lists must be of the same length.")

            err : list = []
            for i in range(len(trueValue)):
                try:
                    err.append(((trueValue[i] - calcValue[i])**2)/(trueValue[i]**2))
                except ZeroDivisionError:
                    err.append((trueValue[i] - calcValue[i])**2)
        except ZeroDivisionError:
            err : float = (trueValue - calcValue)**2

        return err

    calcValues : dict[str,np.ndarray] = {key: None for k,key in enumerate(list(loads.keys())[3:])}
    strakeErrors : dict[str,float | list] = calcValues.copy()
    check : dict[str,bool | list] = calcValues.copy() # shallow copy is OK

    # P = mean(N_z)*2piR
    calcValues["P"] = (np.amax(N[2,:,-1])+np.amin(N[2,:,-1])) *2*math.pi*r /2
    # Q = 1/2*range(N_zth)*piR
    calcValues["Q"] = (np.amax(N[1,:,-1])-np.amin(N[1,:,-1]))*0.5*math.pi*r
    # T = mean(N_zth)*2pi*R^2
    calcValues["T"] = (np.amax(N[1,:,-1])+np.amin(N[1,:,-1])) *2*math.pi*(r**2) /2
    # M = 1/2*range(N_z)*pi*R^2
    calcValues["M"] = (np.amax(N[2,:,-1])-np.amin(N[2,:,-1]))*0.5 *math.pi*(r**2)

    for k,key in enumerate(list(loads.keys())[3:]):
        loads[key] = float(loads[key]) # change dtype
        strakeErrors[key] = relativeSquaredError(trueValue=loads[key],calcValue=float(np.sum(calcValues[key])))
        check[key] = strakeErrors[key] <= tol

    return check, strakeErrors

# add command line parser
parser = argparse.ArgumentParser('Strake optimiser',
                                    description="Wind Turbine Optimiser",
                                    epilog="Thanks, bye.")
parser.add_argument('-g','--geoms', type=str, required=True, help="Initial strake geometries file, in 'inputs/'")
parser.add_argument('-l','--loads', type=str, required=True, help="Applied loads file, in 'inputs/'")
parser.add_argument('-m','--minimize', type=bool, default=True, help="Boolean: optimize? Default=True")
parser.add_argument('-o','--outfile', type=bool, default=True, help="Boolean: Output geometry and loading. Default=True")
parser.add_argument('-d','--topdiameter', type=float, help="Top of tower's diameter [in metres], otherwise taken from --geoms")

def main() -> None:
    # parse required input arguments
    args = parser.parse_args()
    geoms_filename : str = args.geoms.split("\\")[1]
    loads_filename : str = args.loads.split("\\")[1]

    # constants & parameters
    f_yk : float = 355e6 # characteristic steel strength [N/m2]
    gamma_M1 : float = 1.1 # material partial factor of safety (EN 1993-1-6 Table 4.2)
    E : float = 210e9 # Young's Modulus [N/m2]
    rho : float = 7850 # material density [kg/m3]
    fabClass : str = "A" # fabrication quality class (one of ["A", "B", "C"])
    beta_max : float = math.pi/8 # maximum apex half-angle [rad]

    # import loading
    loadNames : list[str] = ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
    loads : dict[str,float] = {load : importLoadMagnitude(load,loads_filename) for load in loadNames}

    # import geometry
    strakeList = listStrakeIDs(geoms_filename)

    # create strake objects
    strakes : dict[str | int, strake] = {strakeID : strake(strakeID,geoms_filename) for i, strakeID in enumerate(strakeList)}
    theta : np.ndarray[float] = np.linspace(0,2*math.pi,361)

    # to be filled during strake loop
    combinedLoading : dict[dict] = {strakeID : None for strakeID in strakes}
    selfWeight : dict[np.ndarray] = combinedLoading.copy()
    cumulativeTimer : dict = combinedLoading.copy()
    minimumValues : dict[str | int, sc.OptimizeResult] = combinedLoading.copy()

    # loop over strakes for sequential optimisation *if* specified by user
    if args.minimize:
        if args.topdiameter: # if specified by user
            jointRadius : float = args.topdiameter / 2
        else: # take from input file
            args.topdiameter = list(strakes.values())[0].r_top * 2
            jointRadius : float = args.topdiameter / 2
        
        startTime : int = time.perf_counter_ns()

        for strakeID, s in strakes.items():
            # update top radius of strake to enforce geometric continuity from previous strake
            s.update(r_top=jointRadius)

            # save loading applied to this strake
            combinedLoading[strakeID] = loads.copy() # shallow copy of current values is fine

            # coordinate arrays
            z : np.ndarray[float] = np.linspace(0,s.h,100)
            Theta, Z = np.meshgrid(theta,z,indexing="ij")

            # constrained & bounded SLSQP minimisation of volume = fun(t, r_bot)
            minimumValues[strakeID] = sc.minimize(volumeFunction,np.array([s.t,s.r_bot]),args=(s),
                                                method="SLSQP", jac="3-point",
                                                constraints=[{"type":"ineq", "fun" : resistanceCheck, "args" : [s,loads,rho,E,f_yk,fabClass,gamma_M1]},
                                                             {"type":"ineq", "fun" : stressInteractionConstraint, "args" : [s,loads,rho,E,f_yk,fabClass,gamma_M1]}],
                                                bounds=[(0.002,0.1),(s.r_top, s.r_top+(s.h*math.sin(beta_max)))])

            tau_xthRd : float = LS3.findShearBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]
            sigma_xRd : float = LS3.findAxialBucklingStress(E, f_yk, fabClass, gamma_M1, s)[0]

            # add loading from this strake to applied loading for next strake
            selfWeight[strakeID] = p_zStresses(rho*-9.81*s.t,s.h,Z)
            loads["P"] += selfWeight[strakeID][2,0,0]*2*math.pi*s.r # self-weight resultant at z=z0 added to axial load
            loads["M"] += loads["Q"]*s.h # adds moment induced by shear force over this strake's height to baseline moment loading
            jointRadius : float = minimumValues[strakeID].get('x')[1] # saves bottom radius of this strake to enforce geometric continuity in next loop

            cumulativeTimer[strakeID] = time.perf_counter_ns() - startTime

    # write optimiser output to csv files *if* specified by user
    if args.outfile:
        with open(f"results/{loads_filename.split('.')[0]}-{geoms_filename.split('.')[0]}-D1-{args.topdiameter}-new-geometry.csv","x") as fil:
            fil.write(f"ID, Volume [m3], t (mm), d2 (bottom) (mm), d1 (top) (mm), h (mm),\n")
            
            for key, s in strakes.items():
                fil.write(f"{key},")
                fil.write(f"{minimumValues[key].get('fun')},")
                fil.write(f"{abs(minimumValues[key].get('x')[0])*1e3},")
                fil.write(f"{abs(minimumValues[key].get('x')[1])*2e3},")
                fil.write(f"{s.r_top*2e3},")
                fil.write(f"{s.h*1e3},")
                fil.write(f"\n")

        with open(f"results/{loads_filename.split('.')[0]}-{geoms_filename.split('.')[0]}-D1-{args.topdiameter}-new-loading.csv","x") as fil:
            fil.write(f"Strake ID, z0 [m], p_r [Pa], p_th [Pa], p_z [Pa], P [N], Q [N], T [Nm], M [Nm], Self-weight Resultant at z0 [N/m],\n")

            for key, s in strakes.items():
                fil.write(f"{key},")
                fil.write(f"{s.z0},")

                for name in loads:
                    fil.write(f"{combinedLoading[key][name]},")
                
                fil.write(f"{abs(selfWeight[key][2,0,0])},")
                fil.write(f"\n")

        with open(f"results/{loads_filename.split('.')[0]}-{geoms_filename.split('.')[0]}-D1-{args.topdiameter}-optimiser.csv","x") as fil:
            fil.write(f"Strake ID, Success, Message, Iterations, Fun iterations, Jac iterations, Hess iterations, Max constraint violation, Cumulative loop time [ns],\n")

            for key in strakes:
                fil.write(f"{key},")
                fil.write(f"{minimumValues[key].get('success')},")
                fil.write(f"{minimumValues[key].get('message')},")
                fil.write(f"{minimumValues[key].get('nit')},")
                fil.write(f"{minimumValues[key].get('nfev')},")
                fil.write(f"{minimumValues[key].get('njev')},")
                fil.write(f"{minimumValues[key].get('nhev')},")
                fil.write(f"{minimumValues[key].get('maxcv')},")
                fil.write(f"{cumulativeTimer[key]},")
                fil.write(f"\n")

if __name__ == "__main__":
    main()