# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 12:26:43 2024

@author: Cerys Morley
"""
import math
import strake

# finds elastic critical axial buckling stress
def calcElasticCriticalAxialBucklingStress(E: float,C_x: float,r: float,t: float) -> float:
    """
    Uses EN 1993-1-6 D.6 to calculate the elastic critical axial buckling stress.

    E: Young's modulus of the elastic material
    C_x: Critical buckling factor under axial compression (see EN 1993-1-6 Annex D)
    t: Strake thickness
    r: Strake radius
    """
    sigma_xRcr : float = 0.605*E*C_x*t/r
    return sigma_xRcr

# find elastic critical shear buckling stress
def calcElasticCriticalShearBucklingStress(E: float,C_tau: float,omega: float,r: float,t: float) -> float:
    """
    Uses EN 1993-1-6 D.40 to calculate the elastic critical shear buckling stress.

    E: Young's modulus of the elastic material
    C_tau: Critical buckling factor under shear (see EN 1993-1-6 Annex D)
    t: Strake thickness
    r: Strake radius
    omega: Relative length of strake (see EN 1993-1-6 D.1)
    """
    tau_xthRcr : float = 0.75*E*C_tau*math.sqrt(1/omega)*t/r
    return tau_xthRcr

# check if axial buckling even needs to be checked
def initialAxialBucklingCheck(E: float,C_x: float,f_yk: float,r:float,t:float) -> bool:
    """
    Checks whether axial buckling stress needs to be evaluated using EN 1993-1-6 D.9.

    E: Young's modulus of the elastic material
    C_x: Critical buckling factor under axial compression (see EN 1993-1-6 Annex D)
    t: Strake thickness
    r: Strake radius
    f_yk: Characteristic yield stress of material
    """
    check : bool = r/t <= C_x*E/(165*f_yk)
    return check

def initialShearBucklingCheck(E: float,f_yk: float,r: float,t: float) -> bool:
    """
    Checks whether shear buckling stress needs to be evaluated using EN 1993-1-6 D.54.

    E: Young's modulus of the elastic material
    t: Strake thickness
    r: Strake radius
    f_yk: Characteristic yield stress of material
    """
    check : bool = r/t <= 0.17*((E/f_yk)**0.67)
    return check

def calcRelativeSlenderness(f_yk: float,sigma_Rcr: float) -> float:
    """
    Calculates the relative slenderness $\bar{\lambda}$ of the shell using
    EN 1993-1-6 Equations 9.19, 9.20, 9.21.

    f_yk: Characteristic yield stress of material. (If evaluating shear stress, use 'f_yk/sqrt(3)')
    sigma_Rcr: elastic critical buckling stress
    """
    lambda_bar = math.sqrt(f_yk/sigma_Rcr)
    return lambda_bar

def calcBucklingReductionFactor(lambda_bar: float,lambda_bar0: float,chi_h: float,alpha: float,beta: float,eta: float) -> float:
    """
    Calculates the elastic-plastic buckling reduction factor using
    EN 1993-1-6 Equations 9.22, 9.23, 9.24.

    lambda_bar: Relative slenderness
    lambda_bar0: 
    chi_h: 
    alpha:
    beta:
    eta:
    """
    lambda_barp : float = calcLambdaBarP(alpha, beta)

    if lambda_bar <= lambda_bar0:
        chi : float = chi_h - (lambda_bar/lambda_bar0)*(chi_h-1)
    elif lambda_bar < lambda_barp:
        chi : float = 1 - beta*(((lambda_bar-lambda_bar0)/(lambda_barp-lambda_bar0))**eta)
    else:
        chi : float = alpha/(lambda_bar**2)

    return chi

def calcImperfectionAmplitude(Q: float,r: float,t: float) -> float:
    """
    Calculates imperfection amplitude $\delta/t$, using EN 1993-1-6 Equations D.14, D.49.
    Works for both axial compression and shear.

    Q: Fabrication quality parameter
    t: Strake thickness
    r: Strake radius
    """
    dt : float = math.sqrt(r/t)/Q
    return dt

def calcAlpha_x(dt: float) -> float:
    """
    Calculates the axial elastic imperfection reduction factor $\alpha_x$,
    using EN 1993-1-6 Equations D.11, D.12, D.13.

    dt: Imperfection amplitude $\delta/t$
    """
    alpha_x = 0.83/(1 + 2.2*(dt**0.75))
    return alpha_x

# calculate shear elastic imperfection reduction factor
def calcAlpha_tau(dt: float):
    # EN 1993-1-6 D.47, D.48
    alpha_tau = 0.96/(1 + 0.5*dt)
    return alpha_tau

# calculate axial plastic range factor
def calcBeta_x(dt: float):
    # EN 1993-1-6 D.15
    beta_x = 1 - (0.75/(1 + 1.1*dt))
    return beta_x

# plastic limit relative slenderness
def calcLambdaBarP(alpha: float,beta: float):
    LP = math.sqrt(alpha/(1-beta)) # EN 1993-1-6 9.25
    return LP

# calculate interaction exponent eta_x0
def calcEta_x0(dt: float):
    eta_x0 = 1.35 - 0.1*dt # EN 1993-1-6 D.16
    return eta_x0

# calculate plastic (?) interaction exponent eta_xp
def calcEta_xp(dt: float):
    eta_xp = 1/(0.45 + 0.72*dt) # EN 1993-1-6 D.17
    return eta_xp

# interaction exponent
def calcEta(lambda_bar: float,lambda_bar0: float,lambda_barp: float,eta_0: float,eta_p: float):
    # EN 1993-1-6 9.26
    eta = ( (lambda_bar*(eta_p - eta_0)) + (lambda_barp*eta_0) - (lambda_bar0*eta_p) )/(lambda_barp-lambda_bar0)
    return eta

# characteristic & design buckling stresses
def calcDesignBucklingStress(chi: float,f_yk: float,gamma_M1: float):
    """Calculates design resistance from material strength."""
    # if doing shear stress, input f_yk = f_yk/sqrt(3)
    sigma_Rk = chi*f_yk # EN 1993-1-6 9.27, 9.28, 9.29

    sigma_Rd = sigma_Rk/gamma_M1 # EN 1993-1-6 9.30, 9.31, 9.32

    return sigma_Rd

# check interactions between stress components
def checkStressInteractions(sigma_Ed: list[float],sigma_Rd: list[float],chi: list[float]) -> bool:
    """Returns binary stress interaction check. Input lists in [th,xth,x] order."""
    k_i, alpha_i = calck_iAndalpha_i(chi)
    # EN 1993-1-6 9.5.3 (4)
    if sigma_Ed[2] > 0:
        # if sigma_x,Ed is tensile
        sigma_Ed[2] = 0
    if sigma_Ed[0] < 0:
        # if sigma_th,Ed is tensile
        sigma_Ed[0] = 0

    # EN 1993-1-6 9.36
    check : bool = ((sigma_Ed[2]/sigma_Rd[2])**k_i[2]) - alpha_i*(sigma_Ed[2]/sigma_Rd[2])*(sigma_Ed[0]/sigma_Rd[0]) + ((sigma_Ed[0]/sigma_Rd[0])**k_i[0]) + ((sigma_Ed[1]/sigma_Rd[1])**k_i[1]) <= 1.0
    return check

# buckling interaction parameters
def calck_iAndalpha_i(chi: list[float]) -> tuple[list[float],float]:
    """Uses chi as list of size 3, [th,xth,x] order to calculate k_i and alpha_i."""
    k_i = []
    k_i.append(1.0 + chi[0]**2) # EN 1993-1-6 9.37
    k_i.append(1.5 + 0.5*(chi[1]**2)) # EN 1993-1-6 9.39
    k_i.append(1.0 + chi[2]**2) # EN 1993-1-6 9.38

    alpha_i = (chi[2]*chi[0])**2 # EN 1993-1-6 9.40

    return k_i, alpha_i

def findAxialBucklingStress(E: float,f_yk: float, fabClass: str, gamma_M1: float,s: strake.strake):
    # axial parameters
    Q_x = {"A":40, "B":25, "C": 16} # EN 1993-1-6 Table D.1
    Q_x_val = Q_x[fabClass]
    chi_xh = 1.10 # EN 1993-1-6 D.19
    lambda_x0 = 0.10 # EN 1993-1-6 D.10
    # check = initialAxialBucklingCheck(E, C_x, f_yk, s.r, s.t)
    # if check == True:
    #     print("Axial buckling does not need to be checked.")

    sigma_xRcr = calcElasticCriticalAxialBucklingStress(E, s.C_x,s.r,s.t)

    lambda_bar = calcRelativeSlenderness(f_yk, sigma_xRcr)

    dt = calcImperfectionAmplitude(Q_x_val, s.r, s.t)
    alpha = calcAlpha_x(dt)
    beta = calcBeta_x(dt)

    eta_x0 = calcEta_x0(dt)
    eta_xp = calcEta_xp(dt)
    lambda_barp = calcLambdaBarP(alpha, beta)

    eta = calcEta(lambda_bar, lambda_x0, lambda_barp, eta_x0, eta_xp)
    chi = calcBucklingReductionFactor(lambda_bar, lambda_x0, chi_xh, alpha, beta, eta)

    sigma_xRd = calcDesignBucklingStress(chi, f_yk, gamma_M1)

    return sigma_xRd, chi

def findShearBucklingStress(E: float, f_yk: float, fabClass: str, gamma_M1:float, s: strake.strake):
    # shear parameters
    Q_tau = {"A":40, "B":25, "C": 16} # EN 1993-1-6 Table D.7
    Q_tau_val = Q_tau[fabClass]
    chi_tauh = 1.0 # EN 1993-1-6 D.53
    lambda_tau0 = 0.40 # EN 1993-1-6 D.50
    beta_tau = 0.60 # EN 1993-1-6 D.51
    eta_tau = 1.0 # EN 1993-1-6 D.52
    # check = initialShearBucklingCheck(E, f_yk, s.r, s.t)
    # if check == True:
    #     print("Axial buckling does not need to be checked.")

    tau_xthRcr = calcElasticCriticalShearBucklingStress(E, s.C_tau, s.omega, s.r, s.t)
    lambda_bar = calcRelativeSlenderness(f_yk/math.sqrt(3), tau_xthRcr)

    dt = calcImperfectionAmplitude(Q_tau_val, s.r, s.t)
    alpha_tau = calcAlpha_tau(dt)

    chi = calcBucklingReductionFactor(lambda_bar, lambda_tau0, chi_tauh, alpha_tau, beta_tau, eta_tau)

    tau_xthRd = calcDesignBucklingStress(chi, f_yk/math.sqrt(3), gamma_M1)

    return tau_xthRd, chi

def findCircumferentialBucklingStress(E: float, f_yk: float, fabClass: str, gamma_M1:float, s: strake.strake) -> tuple[float,float]:
    return 1.0, 0.0
    # raise NotImplementedError("Functionality not implemented.")

# if running this file directly
if __name__ == "__main__":
    # check that strake doesn't buckle
    # constants
    f_yk = 345e6 # characteristic steel strength [N/m2]
    gamma_M1 = 1.1 # EN 1993-1-6 Table 4.2
    E = 210e9 # Young's Modulus [N/m2]

    # geometry import
    filename = "Sadowskietal2023-benchmarkgeometries.csv"
    strakeID = 107
    s = strake.strake(strakeID,filename)

    # find design resistances
    sigma_xRd, chi_x = findAxialBucklingStress(E,f_yk,"A",gamma_M1,s)
    tau_xthRd, chi_tau = findShearBucklingStress(E,f_yk,"A",gamma_M1,s)