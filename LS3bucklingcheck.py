# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 12:26:43 2024

@author: Cerys Morley
"""
import math
import conicalgeometrytocylindrical as geometry

# finds elastic critical axial buckling stress
def calcElasticCriticalAxialBucklingStress(E,C_x,t,r):
    # EN 1993-1-6 D.6
    sigma = 0.605*E*C_x*t/r;
    return sigma

# find critical buckling factor C_x
def calcC_x(omega,r,t):

    if omega < 1.7: # EN 1993-1-6 D.3
        # short length
        C_x = 1.36 - (1.83/omega) + (2.07/(omega**2)); # EN 1993-1-6 D.8

    elif omega < 1.43*r/t: # EN 1993-1-6 D.4
        # medium length
        C_x = 1; # EN 1993-1-6 D.7

    else:
        # long length
        C_x = 1;

    return C_x

# relative length of shell segment
def calcRelativeLength(L,r,t):
    # EN 1993-1-6 Equation D.1
    omega = L/math.sqrt(r*t);
    return omega

# relative slenderness of shell segment
def calcRelativeSlenderness(f_yk,sigma):
    # EN 1993-1-6 Equations 9.19, 9.20, 9.21
    # list in x, theta, xtheta order
    # sigma must be in the same order & of size (3,)
    f_y = [f_yk, f_yk, f_yk/math.sqrt(3)];

    lambda_bar = math.sqrt(f_y/sigma);

    return lambda_bar

# elastic-plastic buckling reduction factor
def calcBucklingReductionFactor(lambda_bar,lambda_bar0,chi_h,alpha,beta,eta):
    lambda_barp = calcLambdaBarP(alpha, beta);

    if lambda_bar <= lambda_bar0:
        chi = chi_h - (lambda_bar/lambda_bar0)*(chi_h-1); # EN 1993-1-6 9.22
    elif lambda_bar < lambda_barp:
        chi = 1 - beta*(((lambda_bar-lambda_bar0)/(lambda_barp-lambda_bar0))**eta); # EN 1993-1-6 9.23
    else:
        chi = alpha/(lambda_bar**2); # EN 1993-1-6 9.24
    return chi

# plastic limit relative slenderness
def calcLambdaBarP(alpha,beta):
    LP = math.sqrt(alpha/(1-beta)); # EN 1993-1-6 9.25
    return LP

# interaction exponent
def calcEta(lambda_bar,lambda_bar0,lambda_barp,eta_0,eta_p):
    # EN 1993-1-6 9.26
    eta = ( (lambda_bar*(eta_p - eta_0)) + (lambda_barp*eta_0) - (lambda_bar0*eta_p) )/(lambda_barp-lambda_bar0);
    return eta

# characteristic & design buckling stresses
def calcDesignBucklingStress(chi,f_yk,gamma_M1):

    f_y = [f_yk, f_yk, f_yk/math.sqrt(3)];

    sigma_Rk = chi*f_y; # EN 1993-1-6 9.27, 9.28, 9.29

    sigma_Rd = sigma_Rk/gamma_M1; # EN 1993-1-6 9.30, 9.31, 9.32
    return sigma_Rd

# check individual components do not exceed design values
def checkIndividualStresses(sigma_Ed,sigma_Rd):
    # gives binary true/false:
    check = sigma_Ed <= sigma_Rd # EN 1993-1-6 9.33, 9.34, 9.35
    return check

# check interactions between stress components
def checkStressInteractions(sigma_Ed,sigma_Rd,k_i,alpha_i):
    # EN 1993-1-6 9.5.3 (4)
    if sigma_Ed[0] < 0:
        # if sigma_x,Ed is tensile
        sigma_Ed[0] = 0;
    if sigma_Ed[1] < 0:
        # if sigma_th,Ed is tensile
        sigma_Ed[1] = 0;

    # gives binary True/False scalar
    # EN 1993-1-6 9.36
    check = ((sigma_Ed[0]/sigma_Rd[0])**k_i[0]) - alpha_i*(sigma_Ed[0]/sigma_Rd[0])*(sigma_Ed[1]/sigma_Rd[1]) + ((sigma_Ed[1]/sigma_Rd[1])**k_i[1]) + ((sigma_Ed[2]/sigma_Rd[2])**k_i[2]) <= 1.0;
    return check

# buckling interaction parameters
def calck_iAndalpha_i(chi):
    # chi is list of size (3,); x, th, xth order
    k_i = [];
    k_i.append(1.0 + chi[0]**2); # EN 1993-1-6 9.37
    k_i.append(1.0 + chi[1]**2); # EN 1993-1-6 9.38
    k_i.append(1.5 + 0.5*(chi[2]**2)); # EN 1993-1-6 9.39

    alpha_i = (chi[0]*chi[1])**2; # EN 1993-1-6 9.40

    return k_i, alpha_i

# check that strake doesn't buckle
filename = "Sadowskietal2023-benchmarkgeometries.csv";
strakeID = 112;

f_yk = 275; # characteristic steel strength [N/mm2]
E = 210e9; # Young's Modulus [N/mm2]

L, r, t = geometry.findStrakeGeometry(filename, strakeID);

omega = calcRelativeLength(L, r, t);
C_x = calcC_x(omega, r, t);

# add if under axial compression statement here?
sigma_xRcr = calcElasticCriticalAxialBucklingStress(E, C_x, t, r);
sigma_thRcr = 0; # expand? not necessary atm.
tau_xthRcr = 0; # expand?
sigma = [sigma_xRcr, sigma_thRcr, tau_xthRcr];

lambda_bar = calcRelativeSlenderness(f_yk, sigma);


lambda_barp = calcLambdaBarP(alpha, beta);
eta = calcEta(lambda_bar, lambda_bar0, lambda_barp, eta_0, eta_p);
chi = calcBucklingReductionFactor(lambda_bar, lambda_bar0, chi_h, alpha, beta, eta);
