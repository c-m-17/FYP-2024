"""
Gradient descent function.

@author: Cerys Morley
"""
import numpy as np
import math
from scipy import optimize as sc

import strake
import LS3bucklingcheck as LS3
import fundamentalstressresultants as fsr

def gradientDescent(x, f, grad, eta: float, alpha: float, i_max: int, tol: float):
    """Steepest gradient descent algorithm, returns iteration history"""

    # initial parameters
    i = 0
    df = 1e10
    dx = 0

    xLog = x
    fLog = f(x)

    while i < i_max and df > tol:

        dx = -eta*grad(x) + alpha*dx
        x += dx

        xLog = np.vstack((xLog,x))
        fLog = np.vstack((fLog,f(x)))

        i += 1
        df = np.absolute(fLog[-1] - fLog[-2])

    return xLog, fLog

def objectiveFunction(t: list, s: strake.strake, rho: float, loads: dict, i: int):
    # set up
    theta = np.linspace(0,2*math.pi,361)
    z = np.linspace(0,s.h,100)
    Theta, Z = np.meshgrid(theta,z,indexing="ij")
 
    selfWeight = np.array(fsr.p_zStresses(rho*-9.81*s.t,s.h,Z)) # selfweight as p_z

    # stresses
    N = fsr.cumulativeStresses(loads,Z,Theta,s)
    N += selfWeight

    sigma_Ed = np.amin(N[i,:,:])/s.t

    # check which component is being checked
    if i == 1:
        sigma_Rd = LS3.findShearBucklingStress(E, f_yk, Q_tau, lambda_tau0, chi_tauh, beta_tau, eta_tau,gamma_M1,s)[0]
    elif i == 2:
        sigma_Rd = LS3.findAxialBucklingStress(E, f_yk, Q_x, lambda_x0, chi_xh,gamma_M1,s)[0]

    f = LS3.checkIndividualStresses(sigma_Ed,sigma_Rd)[1]

    return f

if __name__ == "__main__":
    
    # objective function
    def f(x):
        return np.sum(x*x + 2)
    
    def grad(x):
        return 2*x
    
    eta = 0.05
    alpha = 0.8
    i_max = 30

    x = np.array([10.1,10.1])

    xLog, fLog = gradientDescent(x,f,grad,eta,alpha,i_max, 0.01)

    print(xLog)
    print(fLog)






