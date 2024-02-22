# -*- coding: utf-8 -*-
"""
This script finds the membrane stress resultants N acting on a cylindrical shell.
The geometry, applied loading and material properties are taken as inputs.
Only *membrane* equilibrium is considered.

@author: Cerys Morley
Created on 19 Feb 2024
"""
## import packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

import conicalgeometrytocylindrical as geometry

## function definitions
### dN over general 2D element with dA=R*dth*dz
def p_rStresses(p_r,R,th_n,z_n):
    dN_th = -p_r*R;
    dN_zth = np.zeros((th_n,z_n));
    dN_z = np.zeros((th_n,z_n));
    return [dN_th, dN_zth, dN_z]

def p_thStresses(p_th,z,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = -np.trapz(p_th,z,axis=1);
    dN_zth = np.meshgrid(z,dN_zth)[1];
    dN_z = np.zeros((th_n,z_n));
    return [dN_th, dN_zth, dN_z]

def p_zStresses(p_z,z,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = np.zeros((th_n,z_n));
    dN_z = -np.trapz(p_z,z,axis=1);
    dN_z = np.meshgrid(z,dN_z)[1];
    return [dN_th, dN_zth, dN_z]

def PStresses(P,R,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = np.zeros((th_n,z_n));
    dN_z = P*np.ones((th_n,z_n))/(2*np.pi*R);
    return [dN_th, dN_zth, dN_z]

def TStresses(T,R,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = T*np.ones((th_n,z_n))/(2*np.pi*R*R);
    dN_z = np.zeros((th_n,z_n));
    return [dN_th, dN_zth, dN_z]

def MStresses(M,R,coords,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = np.zeros((th_n,z_n));
    dN_z = (M/(np.pi*R*R))*np.cos(coords[1]);
    return [dN_th, dN_zth, dN_z]

def QStresses(Q,R,coords,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = Q*np.sin(coords[1])/(np.pi*R);
    dN_z = np.zeros((th_n,z_n));
    return [dN_th, dN_zth, dN_z]

### from other Ns
def N_thInducesN_zth(N_th,R,th,z,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = -np.trapz(np.gradient(N_th,th,axis=0)/R,z,axis=1);
    dN_zth = np.meshgrid(z,dN_zth)[1];
    dN_z = np.zeros((th_n,z_n));
    return [dN_th, dN_zth, dN_z]

def N_zthInducesN_z(N_zth,R,th,z,th_n,z_n):
    dN_th = np.zeros((th_n,z_n));
    dN_zth = np.zeros((th_n,z_n));
    dN_z = -np.trapz(np.gradient(N_zth,th,axis=0)/R,z,axis=1);
    dN_z = np.meshgrid(z,dN_z)[1];
    return [dN_th, dN_zth, dN_z]


## input cylindrical shell geometry
# - radius "R" (mm)
# - height "H" (mm)
# - thickness "t" (mm)
# - segment angle "TH" (radians)
TH = 2*np.pi;

filename = "Sadowskietal2023-benchmarkgeometries.csv";
strakeID = int(104);
# gives list of [H,R,t] values for given strake, in mm
geom = geometry.findStrakeGeometry(filename, strakeID);
H = geom[0][0];
R = geom[1][0];
t = geom[2][0];

## cylindrical coordinates r,th,z (radial, circumferential, meridional)
### arrays where th = row, z = column
th_n = 360;
z_n = math.ceil(H/2);
th = np.linspace(0,TH,th_n);
z = np.linspace(0,H,z_n);
### coords[0] = z-coords, coords[1] = th-coords
coords = np.meshgrid(z,th);


## input applied loading
### loads = f(th,z) *only*, no variation in r considered.

#### - pressures p_r, p_th, p_z (N/mm2 (= MPa))
p_r = (1 + 1*(1-np.divide(coords[0],H)));
p_th = 0.1*np.ones((th_n,z_n));
p_z = -0.5*p_r;

#### - axial force P (N)
P = -1000; # (+ve z-dir)

#### - uniform torque T (Nmm (= kNm))
T = 1; # (+ve th-dir)

#### - global bending moment M (Nmm (=kNm))
M = 2; # (C around y-axis)

#### - transverse shear force Q (N)
Q = 5; # (+ve x-dir)


## find global loading
# reaction forces R(th,z=0) [N/mm]
R_z = -P*np.ones((th_n,))/(2*R*np.pi); # from axial force
R_z += -np.trapz(p_z,z,axis=1); # from tangential pressure

R_th = -Q*np.sin(th)/(2*np.pi*R); # from shear force Q
R_th += -T*np.ones((th_n,))/(2*np.pi*R*R); # from torsion
R_th += -np.trapz(p_th,z,axis=1); # from circ pressure

R_rho = -Q*np.cos(th)/(2*np.pi*R); # from shear force Q
R_rho += -np.trapz(p_r,z,axis=1); # from normal pressure

R_M = -(M+ Q*H); # Nmm (C)

# internal force distributions f(th,z)
SFD = np.meshgrid(z,-R_rho)[1]; # shear force (radial) (N)
AFD = np.meshgrid(z,-R_z)[1]; # axial force (N)
BMD = -R_M-(Q*coords[0]); # bending moment (Nmm)
GTD = np.meshgrid(z,-R_th)[1]; # torque diagram (Nmm)

## plotted
plt.figure(1);
plt.plot(AFD[0,:],z);

plt.figure(2);
plt.plot(SFD[0,:],z);

plt.figure(3);
plt.plot(BMD[0,:],z);

plt.figure(4);
plt.plot(GTD[0,:],z);



## superimpose Ns
### from applied loading, N[0] = N_th, N[1] = N_zth, N[2] = N_z
N = np.array([np.empty((th_n,z_n))]*3); # initialise
N += np.asarray(p_rStresses(p_r, R,th_n,z_n),dtype=float);
N += np.asarray(p_thStresses(p_th,z,th_n,z_n),dtype=float);
N += np.asarray(p_zStresses(p_z,z,th_n,z_n),dtype=float);
N += np.asarray(PStresses(P,R,th_n,z_n),dtype=float);
N += np.asarray(TStresses(T,R,th_n,z_n),dtype=float);
N += np.asarray(MStresses(M,R,coords,th_n,z_n),dtype=float);
N += np.asarray(QStresses(Q,R,coords,th_n,z_n),dtype=float);

### dependency on other Ns
N += np.asarray(N_thInducesN_zth(N[0,:,:],R,th,z,th_n,z_n),dtype=float);

#### add boundary function f_1 first before assessing N_zth contribution to N_z
f_1 = -N[1,:,0]; # N_zth(z=0) + f_1 = 0
N[1,:,:] += np.meshgrid(z,f_1)[1];

N += np.asarray(N_zthInducesN_z(N[1,:,:], R, th, z, th_n, z_n),dtype=float);

#### add boundary function f_2
f_2 = -N[2,:,-1]; # N_z(z=H) = 0 (free edge) 
N[2,:,:] += np.meshgrid(z,f_2)[1];

## plot global loading
# - use matplotlib
# - function to add standard plot formatting
# - plot AFD on z-axis
plt.plot(z,N[0,0,:]);
plt.plot(th,N[1,:,0]);
plt.plot(z,N[2,0,:]);

# - plot SFD on z-axis
# - plot BMD on z-axis
# - plot torque on r-t plane?

# ## membrane equilibrium
# # - radial equilibrium
# N_th = p_r*R;

# # - circumferential equilbrium
# dp_rdth = np.gradient(p_r,axis=0);
# N_zth = -np.trapz(dp_rdth+p_th,axis=1);
# N_zth = np.transpose(N_zth*np.ones((z_n,th_n)));
# f_1 = -N_zth; # N_zt(z=H) = 0 (BC1)
# N_zth += f_1;

# # - meridional equilibrium
# dN_zthdth = np.gradient(N_zth,axis=0);
# N_z = -np.trapz(dN_zthdth/R + p_z,axis=1)
# f_2 = -N_z; # N_z(z=H) = 0 (BC2)
# N_z += f_2;

# plot N distribution
# - plot on z-axis

# find critical N
# - max()

# identify failure mechanism
# - bursting or buckling?
# - location
