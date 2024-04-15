"""
Error checking functions
@author: Cerys Morley
"""
import numpy as np
import math

def absoluteRelativeError(trueValue, calcValue):
    """Error value."""
    if trueValue != 0:
        err = abs(trueValue - calcValue)/abs(trueValue)
    else:
        err = 0
    return err

# functions
def printBoolCheck(f, check: str,binary: bool,val: float):

    if binary:
        f.write(check +" passed.\n")
    else:
        f.write(check +" failed.\n")

    f.write("Value: "+str(val)+"\n")

def calculateStressErrors(N: np.ndarray[float], loads: dict[str,float], r: float, tol: float):
    """Calculates error between numerical and analytical membrane stresses to return binary tolerance check"""

    calcValues = {key: [] for k,key in enumerate(list(loads.keys())[3:])}
    strakeErrors = {key: None for k,key in enumerate(list(loads.keys())[3:])}
    check = strakeErrors.copy()

    # P = mean(N_z)*2piR
    calcValues["P"].append((np.amax(N[2,:,-1])+np.amin(N[2,:,-1])) *2*math.pi*r /2)
    # Q = 1/2*range(N_zth)*piR
    calcValues["Q"].append((np.amax(N[1,:,-1])-np.amin(N[1,:,-1]))*0.5*math.pi*r)
    # T = mean(N_zth)*2pi*R^2
    calcValues["T"].append((np.amax(N[1,:,-1])+np.amin(N[1,:,-1])) *2*math.pi*(r**2) /2)
    # M = 1/2*range(N_z)*pi*R^2
    calcValues["M"].append((np.amax(N[2,:,-1])-np.amin(N[2,:,-1]))*0.5 *math.pi*(r**2))

    for k,key in enumerate(list(loads.keys())[3:]):
        strakeErrors[key] = absoluteRelativeError(loads[key],sum(calcValues[key]))
        check[key] = strakeErrors[key] <= tol

    return [check, strakeErrors]