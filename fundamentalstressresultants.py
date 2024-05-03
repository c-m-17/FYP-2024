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

def cumulativeStresses(loads: dict[str,float], Z: np.ndarray[float,float], Theta: np.ndarray[float,float], s: strake.strake) -> np.ndarray:
    """
    Sums the contributions to membrane stress resultants
    from all load types.

    loads: constant values of different applied load types named as exactly [p_r, p_th, p_z, P, Q, T, M]
    Z: Local z-coordinate values corresponding to strake, in interval [0,strake height]
    Theta: Local th-coordinate values corresponding to strake, in interval [0,2*pi]
    s: strake object
    
    returns: cumulative N_theta, N_ztheta, N_z distributions as 3D numpy array where indices
    correspond to [N component, theta-coord, z-coord].
    """

    N = np.zeros((3,*np.shape(Z)))

    # add stresses
    N += p_rStresses(loads["p_r"],s.r,Z)
    N += p_thStresses(loads["p_th"],Z)
    N += p_zStresses(loads["p_z"],s.h,Z)
    N += PStresses(loads["P"],s.r,Z)
    N += QStresses(loads["Q"],s.r,s.h-Z,Theta)
    N += TStresses(loads["T"],s.r,Theta)
    N += MStresses(loads["M"],s.r,Theta)

    return N