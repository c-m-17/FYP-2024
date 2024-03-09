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
import math
import matplotlib.pyplot as plt
from matplotlib import cm

import conicalgeometrytocylindrical as geometry
import globalloads as gl
from plotter import plotter

## function definitions
### dN over general 2D element with dA=R*dth*dz
def p_rStresses(p_r,R,Z):
    dN_th = p_r*R*np.ones(np.shape(Z));
    dN_zth = np.zeros(np.shape(Z));
    dN_z = np.zeros(np.shape(Z));
    return np.array([dN_th, dN_zth, dN_z])

def p_thStresses(p_th,Z):
    dN_th = np.zeros(np.shape(Z));
    dN_zth = -p_th*Z;
    dN_z = np.zeros(np.shape(Z));
    return np.array([dN_th, dN_zth, dN_z])

def p_zStresses(p_z,Z):
    dN_th = np.zeros(np.shape(Z));
    dN_zth = np.zeros(np.shape(Z));
    dN_z = -p_z*Z;
    return np.array([dN_th, dN_zth, dN_z])

def PStresses(P,R,Z):
    dN_th = np.zeros(np.shape(Z));
    dN_zth = np.zeros(np.shape(Z));
    dN_z = P*np.ones(np.shape(Z))/(2*math.pi*R);
    return np.array([dN_th, dN_zth, dN_z])

def TStresses(T,R,Theta):
    dN_th = np.zeros(np.shape(Theta));
    dN_zth = T*np.ones(np.shape(Theta))/(2*math.pi*R**2);
    dN_z = np.zeros(np.shape(Theta));
    return np.array([dN_th, dN_zth, dN_z])

def MStresses(M,R,Theta):
    dN_th = np.zeros(np.shape(Theta));
    dN_zth = np.zeros(np.shape(Theta)); # *applied* moment gradient = 0. gradient from Q included there.
    dN_z = M*np.cos(Theta)/(math.pi*R**2); # cos(th)
    return np.array([dN_th, dN_zth, dN_z])

def QStresses(Q,R,Theta,lever_arm):
    dN_th = np.zeros(np.shape(Theta));
    dN_zth = Q*np.sin(Theta)/(math.pi*R); # sin(th)
    dN_z = -Q*np.cos(Theta)*lever_arm/(math.pi*R**2); # cos(th)
    return np.array([dN_th, dN_zth, dN_z])

## input cylindrical shell geometry
# - radius "R" (mm)
# - strake height "h" (mm)
# - thickness "t" (mm)

geoms_filename = "Sadowskietal2023-benchmarkgeometries.csv";
strakeID = 114;

# gives list of [H,R,t] values for given strake, in mm
h, R, t = geometry.findStrakeGeometry(geoms_filename, strakeID);

z0, H = geometry.findStrakePositionGlobal(geoms_filename, strakeID);

## input global loading
loads_filename = "Sadowskietal2023-benchmarkloads-LC1.csv";

### pressure loads (N/mm2 = (MPa))
p_r = gl.importLoadMagnitude("p_r",loads_filename);
p_th = gl.importLoadMagnitude("p_th",loads_filename);
p_z = gl.importLoadMagnitude("p_z",loads_filename);

#### - axial force P (N)
P = gl.importLoadMagnitude("P",loads_filename); # (+ve z-dir)
#### - transverse shear force Q (N)
Q = gl.importLoadMagnitude("Q",loads_filename); # (+ve x-dir)
#### - uniform torque T (Nmm (= kNm))
T = gl.importLoadMagnitude("T",loads_filename); # (+ve th-dir)
#### - global bending moment M (Nmm (=kNm))
M = gl.importLoadMagnitude("M",loads_filename); # (C around y-axis)

## cylindrical coordinates r,th,z (radial, circumferential, meridional)
th_n = 361; sec_angle = 2*math.pi;
Theta, Z = np.meshgrid(np.linspace(0,sec_angle,th_n), np.linspace(z0,h,100),indexing="ij");
X, Y = R*np.cos(Theta), R*np.sin(Theta); # for plotting only

## find global internal forces
### reaction forces R(th,z*=0) [N/mm]
R_z = -(P/(sec_angle*R) + p_z*(H-z0))*np.ones(np.shape(Theta)); # vertical reaction force
# R_x = - (Q - p_r*(H-z0)*np.cos(theta)); # p_r cancels over 2pi integral
R_th = - (Q*np.sin(Theta)/(sec_angle*R) + p_th*(H-z0));
R_r = - (Q*-np.cos(Theta)/(sec_angle*R) + p_r*(H-z0)) ;

R_T = -(T + p_th*sec_angle*R*R*(H-z0))*np.ones(np.shape(Theta)); # reaction torque [Nmm]
R_M = -(M + Q*(H-z0))*np.ones(np.shape(Theta)); # reaction moment [Nmm]

### internal forces at z*=h
titles = ["Shear (radial) Force","Shear (theta) Force","Axial Force Diagram","Bending Moment Diagram","Global Torque Diagram"];

#### Nz0 i.e. N(z*=0) = -ve reaction forces
Nz0 = np.array([-R_r,-R_th,-R_z,-R_M,-R_T]); # same order as "titles"
diagrams = {k : np.zeros(np.shape(Theta)) for k in titles};

for k,key in enumerate(titles):
    diagrams[key] = Nz0[k,:,:]; # set z*=0 values correctly
    # also initialises other values to be edited below

#### set z*=(0,h] values correctly
diagrams[titles[0]] += -p_r*(Z-z0); # radial shear [N/mm]
diagrams[titles[1]] += -p_th*(Z-z0); # circumferential shear [N/mm]
diagrams[titles[2]] += -p_z*(Z-z0); # axial force (N/mm)
diagrams[titles[3]] += -Q*(Z-z0); # bending moment (Nmm)
diagrams[titles[4]] += -p_th*sec_angle*(R**2)*(Z-z0); # torque (Nmm)

## superimpose Ns for INDIVIDUAL STRAKE
### from applied loading, N[0] = N_th, N[1] = N_zth, N[2] = N_z
N_components = ["$N_\theta$: Circumferential Membrane Stress Resultant","$N_{z\theta}$: Shear Membrane Stress Resultant","$N_z$: Meridional Membrane Stress Resultant"];
N = np.zeros(np.shape(Nz0))[0:3,:,:]; # initialise

# adding fundamental results at all Z=[z0,h]
N += p_rStresses(p_r, R, Z);
N += p_thStresses(p_th, Z);
N += p_zStresses(p_z, Z);
N += PStresses(P, R, Z);
N += TStresses(T, R, Theta);
N += MStresses(M,R,Theta);
N += QStresses(Q,R,Theta,h-Z);

# output values at z*=H
print("\nmax N_th(z*=h) = "+str(max(N[0,:,-1])));
print("\nmax N_zth(z*=h) = "+str(max(N[1,:,-1])));
print("\nmax N_z(z*=h) = "+str(max(N[2,:,-1])));

## error between Ns and diagrams
tol = 1e-3;
print("Global Equilibrium satisfied?\n");
### N_z
absrelerr = abs(P - (np.mean(N[2,:,:])*2*math.pi*R))/abs(P); # base (not th dependent)
print(absrelerr<=tol);
absrelerr = abs(M - (np.amax(N[2,:,-1])-np.amin(N[2,:,-1]))*0.5*math.pi*R*R)/abs(M);
print(absrelerr<=tol);
### N_zth
absrelerr = abs(T - (np.mean(N[1,:,:])*sec_angle*R*R))/abs(T);
print(absrelerr<=tol);
absrelerr = abs(Q - (np.amax(N[1,:,:])-np.amin(N[1,:,:]))*0.5*math.pi*R)/abs(Q);
print(absrelerr<=tol);


# plotting

## plot strake geometry
fig = plt.figure();
ax = fig.add_subplot(projection="3d");
ax.plot_surface(X,Y,Z);
ax.set_aspect('auto');

ax.set_title("Strake "+str(strakeID)+" Geometry");
ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[]);
ax.annotate("Height "+str(h)+" mm,\nRadius "+str(R)+" mm",(50,200),xycoords="figure points");
plt.show();


## plot global loading "diagrams"
for k,key in enumerate(titles):
    fig = plt.figure();
    ax = fig.add_subplot(projection="3d");
    ax.plot_surface(Theta,Z,diagrams[key]);

    ax.set_title(key);
    plt.show();


## plot Ns
for i in range(3):
    fig = plt.figure();
    ax = fig.add_subplot(projection="3d")

    ax.plot_surface(Theta,Z,N[i,:,:],cmap=cm.coolwarm);
    ax.set_title(N_components[i]);
    plt.show();
