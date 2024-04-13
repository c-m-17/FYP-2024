# -*- coding: utf-8 -*-
"""
This script finds the fundamental stress resultants for different types of constant loads.
Units: N, mm

@author: Cerys Morley
2024
"""
import numpy as np
import math
import strake

# dN over general 2D element with dA=R*dth*dz
# see ToS Lecture 2.1 Notes

def p_rStresses(p_r: float, R: float, Z: np.ndarray[float]):
    dN_th = p_r*R*np.ones(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = np.zeros(np.shape(Z))
    return [dN_th, dN_zth, dN_z]

def p_thStresses(p_th: float, Z: np.ndarray[float]):
    dN_th = np.zeros(np.shape(Z))
    dN_zth = -p_th*Z
    dN_z = np.zeros(np.shape(Z))
    return [dN_th, dN_zth, dN_z]

def p_zStresses(p_z: float, H: float, Z: np.ndarray[float]):
    dN_th = np.zeros(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = p_z*(H-Z) # in z-coords
    return [dN_th, dN_zth, dN_z]

def PStresses(P: float, R: float, Z: np.ndarray[float]):
    dN_th = np.zeros(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = P*np.ones(np.shape(Z))/(2*math.pi*R)
    return [dN_th, dN_zth, dN_z]

def TStresses(T: float, R: float, Theta: np.ndarray[float]):
    dN_th = np.zeros(np.shape(Theta))
    dN_zth = T*np.ones(np.shape(Theta))/(2*math.pi*R**2)
    dN_z = np.zeros(np.shape(Theta))
    return [dN_th, dN_zth, dN_z]

def MStresses(M: float, R: float, Theta: np.ndarray[float]):
    dN_th = np.zeros(np.shape(Theta))
    dN_zth = np.zeros(np.shape(Theta)) # *applied* moment gradient = 0. gradient from Q included there.
    dN_z = -M*np.cos(Theta)/(math.pi*R**2) # cos(th)
    return [dN_th, dN_zth, dN_z]

def QStresses(Q: float, R: float, lever_arm: np.ndarray[float], Theta: np.ndarray[float,float]):
    dN_th = np.zeros(np.shape(Theta))
    dN_zth = Q*np.sin(Theta)/(math.pi*R) # sin(th)
    dN_z = -Q*np.cos(Theta)*lever_arm/(math.pi*R*R) # cos(th)
    return [dN_th, dN_zth, dN_z]

def cumulativeStresses(loads: dict[str,float], Z: np.ndarray[float], Theta: np.ndarray[float], s: strake.strake, H: float):
    """Sums the contributions from all load types to the membrane stress resultants"""

    # initialise N: ndarray[float]
    N = np.zeros((3,*np.shape(Z)))

    # add stresses
    N += np.array(p_rStresses(loads["p_r"],s.r,Z))
    N += np.array(p_thStresses(loads["p_th"],Z))
    N += np.array(p_zStresses(loads["p_z"],H,Z))
    N += np.array(PStresses(loads["P"],s.r,Z))
    N += np.array(QStresses(loads["Q"],s.r,H-Z,Theta))
    N += np.array(TStresses(loads["T"],s.r,Theta))
    N += np.array(MStresses(loads["M"],s.r,Theta))

    return N