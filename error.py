"""
Error checking functionality.

@author: Cerys Morley
"""
import numpy as np
import math

def relativeError(*,trueValue, calcValue):
    """
    Calculates the relative error between two values.

    trueValue: exact value, which is also used to normalise error (unless it's equal to zero).
    calcValue: calculated value being checked

    returns: normalised decimal error value
    """
    # if type(trueValue) != type(calcValue):
    #     raise NotImplementedError("Input arguments must have same dtype.")
    try:
        err : float = (trueValue - calcValue)/abs(trueValue)
    except TypeError:
        err : list = []
        for i in range(len(trueValue)):
            err.append((trueValue[i] - calcValue[i])/abs(trueValue[i]))
    except ZeroDivisionError:
        err = (trueValue - calcValue)

    return err

def relativeSquaredError(*,trueValue, calcValue):
    """
    Calculates the relative squared error between two values.

    trueValue: exact value, which is also used to normalise error (unless it's equal to zero).
    calcValue: calculated value being checked

    returns: normalised decimal error value
    """
    # if type(trueValue) != type(calcValue):
    #     raise NotImplementedError("Input arguments must have same dtype.")
    try:
        if trueValue <= 1e-10:
            err = (trueValue - calcValue)**2
        else:
            err : float = ((trueValue - calcValue)**2)/(trueValue**2)
    except TypeError:
        err : list = []
        for i in range(len(trueValue)):
            err.append(((trueValue[i] - calcValue[i])**2)/(trueValue[i]**2))
    except ZeroDivisionError:
        err = (trueValue - calcValue)**2

    return err

def printUtilisationCheck(val: bool | float) -> str:
    """
    Adds 'bool, % rel error' into a string.

    val: value to insert, as either a decimal error or a True/False boolean.
    """
    line : str = ""
    # negative relative error means design < resistance
    if type(val) == bool: # for interaction check
        if val: # success
            line += "True,0,"
        else: # failure
            line += "False,100,"
    elif val >= 0: # success
        line += "True," + str(val*100)+","
    else: # failure
        line += "False," + str(val*100)+","
    
    line += "\n"

    return line

def calculateStressErrors(N: np.ndarray[float], loads: dict[str,float], r: float, tol=1e-10) -> tuple[dict[str,bool],dict[str,float]]:
    """
    Checks global equilibrium by checking theoretical membrane stress resultants against
    calculated ones.

    N: calculated membrane stress resultants in [theta, shear, z] components in that order.
    loads: value of different applied load types in order: [p_r, p_th, p_z, P, Q, T, M]
    r: strake radius
    tol: error tolerance (optional)
    
    returns: Global equilibrium passed? True/False, Relative error value
    """

    calcValues : dict[str,list[float]] = {key: [] for k,key in enumerate(list(loads.keys())[3:])}
    strakeErrors : dict[str,float] = {key: None for k,key in enumerate(list(loads.keys())[3:])}
    check : dict[str,bool] = strakeErrors.copy()

    # P = mean(N_z)*2piR
    calcValues["P"].append((np.amax(N[2,:,-1])+np.amin(N[2,:,-1])) *2*math.pi*r /2)
    # Q = 1/2*range(N_zth)*piR
    calcValues["Q"].append((np.amax(N[1,:,-1])-np.amin(N[1,:,-1]))*0.5*math.pi*r)
    # T = mean(N_zth)*2pi*R^2
    calcValues["T"].append((np.amax(N[1,:,-1])+np.amin(N[1,:,-1])) *2*math.pi*(r**2) /2)
    # M = 1/2*range(N_z)*pi*R^2
    calcValues["M"].append((np.amax(N[2,:,-1])-np.amin(N[2,:,-1]))*0.5 *math.pi*(r**2))

    for k,key in enumerate(list(loads.keys())[3:]):
        strakeErrors[key] = relativeSquaredError(trueValue=float(loads[key]),calcValue=sum(calcValues[key]))
        check[key] = strakeErrors[key] <= tol

    return check, strakeErrors