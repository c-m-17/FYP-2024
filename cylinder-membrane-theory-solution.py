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

import conicalgeometrytocylindrical as geometry
import globalloads as gl
from plotter import plotter

## function definitions
### dN over general 2D element with dA=R*dth*dz
def p_rStresses(p_r,R):
    dN_th = p_r*R;
    dN_zth = 0;
    dN_z = 0;
    return [dN_th, dN_zth, dN_z]

def p_thStresses(p_th,h):
    dN_th = 0;
    dN_zth = -p_th*h;
    dN_z = 0;
    return [dN_th, dN_zth, dN_z]

def p_zStresses(p_z,h):
    dN_th = 0;
    dN_zth = 0;
    dN_z = -p_z*h;
    return [dN_th, dN_zth, dN_z]

def PStresses(P,R):
    dN_th = 0;
    dN_zth = 0;
    dN_z = P/(2*math.pi*R);
    return [dN_th, dN_zth, dN_z]

def TStresses(T,R):
    dN_th = 0;
    dN_zth = T/(2*math.pi*R*R);
    dN_z = 0;
    return [dN_th, dN_zth, dN_z]

def MStresses(M,R,theta):
    # dN_th = 0;
    # dN_zth = 0; # *applied* moment gradient = 0. gradient from Q included there.
    dN_z = M*np.cos(theta)/(math.pi*R*R); # cos(th)
    return dN_z

def QStresses(Q,R,theta,lever_arm):
    # dN_th = 0;
    dN_zth = Q*np.sin(theta)/(math.pi*R); # sin(th)
    dN_z = -Q*np.cos(theta)*lever_arm/(math.pi*R*R); # cos(th)
    return np.array([dN_zth, dN_z])

## input cylindrical shell geometry
# - radius "R" (mm)
# - strake height "h" (mm)
# - thickness "t" (mm)

filename = "Sadowskietal2023-benchmarkgeometries.csv";
strakeID = 104;
# gives list of [H,R,t] values for given strake, in mm
geom = geometry.findStrakeGeometry(filename, strakeID);
h = geom[0];
R = geom[1];
t = geom[2];

z0 = geometry.findStrakePositionGlobal(filename, strakeID)[0];
H = geometry.findStrakePositionGlobal(filename, strakeID)[1];

## input global loading
filename = "Sadowskietal2023-benchmarkloads-LC2.csv";

### pressure loads (N/mm2 = (MPa)) 
p_r = gl.importLoadMagnitude("p_r",filename);
p_th = gl.importLoadMagnitude("p_th",filename);
p_z = gl.importLoadMagnitude("p_z",filename);

#### - axial force P (N)
P = gl.importLoadMagnitude("P",filename); # (+ve z-dir)
#### - transverse shear force Q (N)
Q = gl.importLoadMagnitude("Q",filename); # (+ve x-dir)
#### - uniform torque T (Nmm (= kNm))
T = gl.importLoadMagnitude("T",filename); # (+ve th-dir)
#### - global bending moment M (Nmm (=kNm))
M = gl.importLoadMagnitude("M",filename); # (C around y-axis)

## cylindrical coordinates r,th,z (radial, circumferential, meridional)
th_n = 361; sec_angle = 2*math.pi;
theta = np.linspace(0,sec_angle,th_n);

## find global internal forces
### reaction forces R(th,z*=0) [N/mm]
R_z = -(P/(sec_angle*R) + p_z*(H-z0))*np.ones(len(theta)); # vertical reaction force
# R_x = - (Q - p_r*(H-z0)*np.cos(theta)); # p_r cancels over 2pi integral
R_th = - (Q*np.sin(theta)/(sec_angle*R) + p_th*(H-z0));
R_r = - (Q*-np.cos(theta)/(sec_angle*R) + p_r*(H-z0)) ;

R_T = -(T + p_th*sec_angle*R*R*(H-z0)); # reaction torque [Nmm]
R_M = -(M + Q*(H-z0)); # reaction moment [Nmm]

### internal forces at z*=H
titles = ["Axial Force Diagram","Shear Force Diagram","Bending Moment Diagram","Global Torque Diagram","Shear (theta) Force","Shear (radial) Force"]

# diagrams = {};

# diagrams[titles[0]] = -R_z - (p_z*h); # axial force (N/mm)
# diagrams[titles[1]] = -R_x - (p_r*h*np.cos(theta)); # shear force (N/mm)
# diagrams[titles[2]] = -R_M - Q*h; # bending moment (Nmm)
# diagrams[titles[3]] = -R_T - (p_th*sec_angle*R*R*h); # torque diagram (Nmm)
# diagrams[titles[4]] = -R_th - (p_th*h); # N/mm
# diagrams[titles[5]] = -R_r - (p_r*h); # N/mm

# #### plotted
# decision = input("Do you want to plot figures? Y/N:  ");
# if decision == "Y":
#     # extract the following into separate function file?
#     nrows = 2;
#     ncols = 1;
    
#     xlabel = "z-coordinate [mm]";
#     ylabels = ["[N]",
#                "[N]",
#                "[Nmm]",
#                "[Nmm]"]
    
    # for i,name in enumerate(titles):
    #     plotter(i+1,ncols,nrows,
    #             name,h,xlabel,
    #             diagrams[name],ylabels[i]);
        


## superimpose Ns for INDIVIDUAL STRAKE at z*=h
### from applied loading, N[0] = N_th, N[1] = N_zth, N[2] = N_z
Nh = np.zeros((3,len(theta))); # initialise
# Nh += np.asarray(p_rStresses(p_r,R));
# Nh += np.asarray(p_thStresses(p_th,h));
# Nh += np.asarray(p_zStresses(p_z,h));
Nh += np.meshgrid(np.ones(th_n),np.asarray(PStresses(P,R)))[1];
Nh += np.meshgrid(np.ones(th_n),np.asarray(TStresses(T,R)))[1];
Nh[2,:] += MStresses(M,R,theta);
Nh[range(1,3),:] += QStresses(Q,R,theta,0);

### boundary functions
#### Nz0 i.e. N(z*=0) = reaction forces
Nz0 = np.array([-R_r,-R_th,-R_z]);

#### BC1
f_1 = Nz0[1,:] + p_th*h;
Nh[1,:] += f_1;

#### BC2
f_2 = Nz0[2,:] + np.gradient(Nh[1,:],theta)*h/R + p_z*h;
Nh[2,:] += f_2;

### plot Ns
# if decision == "Y":
#     nrows = 2;
#     ncols = 1;
    
#     ylabel = "[N/mm]";
    
#     titles = ["N_th: Circumferential Membrane Stress Resultant",
#               "N_zth: Shear Membrane Stress Resultant",
#               "N_z: Meridional Membrane Stress Resultant"];
    
#     for i,name in enumerate(titles):
#         plotter(i+5,ncols,nrows,
#                 name,z,xlabel,
#                 N[i]*zshapearray,ylabel);
