# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 12:26:43 2024

@author: Cerys Morley
"""
import math
import strake

# finds elastic critical axial buckling stress
def calcElasticCriticalAxialBucklingStress(E: float,C_x: float,r: float,t: float):
    # EN 1993-1-6 D.6
    sigma_xRcr = 0.605*E*C_x*t/r
    return sigma_xRcr

# find elastic critical shear buckling stress
def calcElasticCriticalShearBucklingStress(E: float,C_tau: float,omega: float,r: float,t: float):
    # EN 1993-1-6 D.40
    tau_xthRcr = 0.75*E*C_tau*math.sqrt(1/omega)*t/r
    return tau_xthRcr

# check if axial buckling even needs to be checked
def initialAxialBucklingCheck(E: float,C_x: float,f_yk: float,r:float,t:float):
    # EN 1993-1-6 D.9
    check = r/t <= C_x*E/(165*f_yk)
    # if true, axial buckling does not need to be checked
    return check

# check if shear buckling evens needs to be checked
def initialShearBucklingCheck(E: float,f_yk: float,r: float,t: float):
    # EN 1993-1-6 D.54
    check = r/t <= 0.17*((E/f_yk)**0.67)
    return check

# relative slenderness of shell segment
def calcRelativeSlenderness(f_yk: float,sigma: float):
    # EN 1993-1-6 9.19, 9.20, 9.21
    lambda_bar = math.sqrt(f_yk/sigma)
    # if doing shear stress, input f_yk = f_yk/sqrt(3)

    return lambda_bar

# elastic-plastic buckling reduction factor
def calcBucklingReductionFactor(lambda_bar: float,lambda_bar0: float,chi_h: float,alpha: float,beta: float,eta: float):
    lambda_barp = calcLambdaBarP(alpha, beta)

    if lambda_bar <= lambda_bar0:
        region = "under squash limit"
        chi = chi_h - (lambda_bar/lambda_bar0)*(chi_h-1) # EN 1993-1-6 9.22
    elif lambda_bar < lambda_barp:
        region = "elastic region"
        chi = 1 - beta*(((lambda_bar-lambda_bar0)/(lambda_barp-lambda_bar0))**eta) # EN 1993-1-6 9.23
    else:
        region = "plastic region"
        chi = alpha/(lambda_bar**2) # EN 1993-1-6 9.24

    return chi, region

# calculate imperfection amplitude
def calcImperfectionAmplitude(Q: float,r: float,t: float):
    # EN 1993-1-6 D.14, D.49
    dt = math.sqrt(r/t)/Q
    # works for axial & shear
    return dt

# calculate axial elastic imperfection reduction factor
def calcAlpha_x(dt: float):
    # EN 1993-1-6 D.11, D.12, D.13
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
    # if doing shear stress, input f_yk = f_yk/sqrt(3)
    sigma_Rk = chi*f_yk # EN 1993-1-6 9.27, 9.28, 9.29

    sigma_Rd = sigma_Rk/gamma_M1 # EN 1993-1-6 9.30, 9.31, 9.32

    return sigma_Rd

# check individual components do not exceed design values
def checkIndividualStresses(sigma_Ed,sigma_Rd):
    "Returns tuple of (bool, error) where error is same size as inputs."
    # gives binary true/false: "Does it exceed buckling stress?"
    check = abs(sigma_Ed) <= abs(sigma_Rd) # EN 1993-1-6 9.33, 9.34, 9.35

    objective = abs(sigma_Rd + sigma_Ed)/abs(sigma_Rd)

    return check, objective

# check interactions between stress components
def checkStressInteractions(sigma_Ed: list[float],sigma_Rd: list[float],chi: list[float]):
    k_i, alpha_i = calck_iAndalpha_i(chi)
    # EN 1993-1-6 9.5.3 (4)
    if sigma_Ed[0] < 0:
        # if sigma_x,Ed is tensile
        sigma_Ed[0] = 0
    if sigma_Ed[1] < 0:
        # if sigma_th,Ed is tensile
        sigma_Ed[1] = 0

    # gives binary True/False scalar
    # EN 1993-1-6 9.36
    check = ((sigma_Ed[0]/sigma_Rd[0])**k_i[0]) - alpha_i*(sigma_Ed[0]/sigma_Rd[0])*(sigma_Ed[1]/sigma_Rd[1]) + ((sigma_Ed[1]/sigma_Rd[1])**k_i[1]) + ((sigma_Ed[2]/sigma_Rd[2])**k_i[2]) <= 1.0
    return [check, 0.0]

# buckling interaction parameters
def calck_iAndalpha_i(chi: list[float]):
    # chi is list of size (3,) x, th, xth order
    k_i = []
    k_i.append(1.0 + chi[0]**2) # EN 1993-1-6 9.37
    k_i.append(1.0 + chi[1]**2) # EN 1993-1-6 9.38
    k_i.append(1.5 + 0.5*(chi[2]**2)) # EN 1993-1-6 9.39

    alpha_i = (chi[0]*chi[1])**2 # EN 1993-1-6 9.40

    return k_i, alpha_i

def findAxialBucklingStress(E: float,f_yk: float, fabClass: str, gamma_M1: float,s: strake.strake):
    # axial parameters
    Q_x = {"Class A":40, "Class B":25, "Class C": 16} # EN 1993-1-6 Table D.1
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
    [chi, region] = calcBucklingReductionFactor(lambda_bar, lambda_x0, chi_xh, alpha, beta, eta)

    sigma_xRd = calcDesignBucklingStress(chi, f_yk, gamma_M1)

    return sigma_xRd, [chi,region]

def findShearBucklingStress(E: float, f_yk: float, fabClass: str, gamma_M1:float, s: strake.strake):
    # shear parameters
    Q_tau = {"Class A":40, "Class B":25, "Class C": 16} # EN 1993-1-6 Table D.7
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

    [chi, region] = calcBucklingReductionFactor(lambda_bar, lambda_tau0, chi_tauh, alpha_tau, beta_tau, eta_tau)

    tau_xthRd = calcDesignBucklingStress(chi, f_yk/math.sqrt(3), gamma_M1)

    return tau_xthRd, [chi,region]

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
    sigma_xRd, chi_x = findAxialBucklingStress(E,f_yk,"Class A",gamma_M1,s)
    tau_xthRd, chi_tau = findShearBucklingStress(E,f_yk,"Class A",gamma_M1,s)