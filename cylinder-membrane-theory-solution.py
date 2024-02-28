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

def p_thStresses(p_th,H):
    dN_th = 0;
    dN_zth = -p_th*H;
    dN_z = 0;
    return [dN_th, dN_zth, dN_z]

def p_zStresses(p_z,H):
    dN_th = 0;
    dN_zth = 0;
    dN_z = -p_z*H;
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

def MStresses(M,R):
    dN_th = 0;
    dN_zth = 0; # *applied* moment gradient = 0. gradient from Q included there.
    dN_z = M/(math.pi*R*R); # cos(th)
    return [dN_th, dN_zth, dN_z]

def QStresses(Q,R,H,z0):
    dN_th = 0;
    dN_zth = Q/(math.pi*R); # sin(th)
    dN_z = -Q*(H-z0)/(math.pi*R*R); # cos(th)
    return [dN_th, dN_zth, dN_z]

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
z_n = 361;
z = np.linspace(0,H,z_n);

## find global force diagrams
### reaction forces R(th,z=0) [N/mm]
R_z = -(P + p_z*H*2*math.pi*R); # vertical reaction force
R_x = - Q; # p_r cancels over 2pi integral
R_T = -(T + p_th*2*math.pi*R*R*H); # reaction torque
R_M = -(M + Q*H); # reaction moment

### internal force distributions f(z), from z=0 to z=H
zshapearray = np.ones((z_n,));
titles = ["Axial Force Diagram","Shear Force Diagram","Bending Moment Diagram","Global Torque Diagram"]

diagrams = dict();

diagrams[titles[0]] = -R_z - (p_z*2*math.pi*R*z); # axial force (N)
diagrams[titles[1]] = -R_x*zshapearray; # shear force (N)
diagrams[titles[2]] = -R_M - Q*z; # bending moment (Nmm)
diagrams[titles[3]] = -R_T - (p_th*2*math.pi*R*R*z); # torque diagram (Nmm)

#### plotted
decision = input("Do you want to plot figures? Y/N:  ");
if decision == "Y":
    # extract the following into separate function file?
    nrows = 2;
    ncols = 1;
    
    xlabel = "z-coordinate [mm]";
    ylabels = ["[N]",
               "[N]",
               "[Nmm]",
               "[Nmm]"]
    
    for i,name in enumerate(titles):
        plotter(i+1,ncols,nrows,
                name,z,xlabel,
                diagrams[name],ylabels[i]);
        


## superimpose Ns for INDIVIDUAL STRAKE
### from applied loading, N[0] = N_th, N[1] = N_zth, N[2] = N_z
N = np.zeros(3); # initialise
N += np.asarray(p_rStresses(p_r, R));
N += np.asarray(p_thStresses(p_th,h));
N += np.asarray(p_zStresses(p_z,h));
N += np.asarray(PStresses(P,R));
N += np.asarray(TStresses(T,R));
N += np.asarray(MStresses(M,R));
N += np.asarray(QStresses(Q,R,H,z0));

### boundary functions
#### BC1
I = list(abs(z - z0)).index(min(abs(z - z0)));

f_1 = (diagrams[titles[1]][I]/(2*math.pi*R)) - N[1];
N[1] += f_1;

# #### BC2
f_2 = (diagrams[titles[0]][I]/(2*math.pi*R)) - N[2];
N[2] += f_2;

# sanity check
internalF = N*2*math.pi*R;
err = [];
for i,name in enumerate(titles):
    if i >= 2:
        break
    err.append(internalF[2-i]/diagrams[name][I]);
    

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
