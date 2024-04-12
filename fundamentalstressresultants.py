# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 13:16:09 2024

@author: Cerys Morley
"""
import numpy as np
import math

# dN over general 2D element with dA=R*dth*dz
# see ToS Lecture 2.1 Notes

def p_rStresses(p_r,R,Z):
    dN_th = p_r*R*np.ones(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = np.zeros(np.shape(Z))
    return np.array([dN_th, dN_zth, dN_z])

def p_thStresses(p_th,Z):
    dN_th = np.zeros(np.shape(Z))
    dN_zth = -p_th*Z
    dN_z = np.zeros(np.shape(Z))
    return np.array([dN_th, dN_zth, dN_z])

def p_zStresses(p_z,Z,H):
    dN_th = np.zeros(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = p_z*(H-Z) # in z-coords
    return np.array([dN_th, dN_zth, dN_z])

def PStresses(P,R,Z):
    dN_th = np.zeros(np.shape(Z))
    dN_zth = np.zeros(np.shape(Z))
    dN_z = P*np.ones(np.shape(Z))/(2*math.pi*R)
    return np.array([dN_th, dN_zth, dN_z])

def TStresses(T,R,Theta):
    dN_th = np.zeros(np.shape(Theta))
    dN_zth = T*np.ones(np.shape(Theta))/(2*math.pi*R**2)
    dN_z = np.zeros(np.shape(Theta))
    return np.array([dN_th, dN_zth, dN_z])

def MStresses(M,R,Theta):
    dN_th = np.zeros(np.shape(Theta))
    dN_zth = np.zeros(np.shape(Theta)) # *applied* moment gradient = 0. gradient from Q included there.
    dN_z = -M*np.cos(Theta)/(math.pi*R**2) # cos(th)
    return np.array([dN_th, dN_zth, dN_z])

def QStresses(Q,R,Theta,lever_arm):
    dN_th = np.zeros(np.shape(Theta))
    dN_zth = Q*np.sin(Theta)/(math.pi*R) # sin(th)
    dN_z = -Q*np.cos(Theta)*lever_arm/(math.pi*R*R) # cos(th)
    return np.array([dN_th, dN_zth, dN_z])
