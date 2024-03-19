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
import fundamentalstressresultants as fsr
import LS3bucklingcheck as LS3

geoms_filename = "Sadowskietal2023-benchmarkgeometries.csv";
loads_filename = "Sadowskietal2023-benchmarkloads-LC1.csv";

## input global loading
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

th_n = 361; sec_angle = 2*math.pi;
# constants
f_yk = 345; # characteristic steel strength [N/mm2]
gamma_M1 = 1.1; # EN 1993-1-6 Table 4.2
E = 210e3; # Young's Modulus [N/mm2]
rho = 7850e-9; # density [kg/mm3]
g = -9.81; # gravity [m/s2]

p_zSW = rho*g; # p_z/t from self-weight [N/mm3]

Q_x = {"Class A":40, "Class B":25, "Class C": 16}; # EN 1993-1-6 Table D.1
chi_xh = 1.10; # EN 1993-1-6 D.19
lambda_x0 = 0.10; # EN 1993-1-6 D.10

Q_tau = {"Class A":40, "Class B":25, "Class C": 16}; # EN 1993-1-6 Table D.7
chi_tauh = 1.0; # EN 1993-1-6 D.53
lambda_tau0 = 0.40; # EN 1993-1-6 D.50
beta_tau = 0.60; # EN 1993-1-6 D.51
eta_tau = 1.0; # EN 1993-1-6 D.52

# find number of strakes
strakeList, H, V = geometry.listStrakeIDs(geoms_filename);

titles = ["Shear (radial) Force","Shear (theta) Force","Axial Force Diagram","Bending Moment Diagram","Global Torque Diagram"];
N_components = ["$N_{\theta}$: Circumferential Membrane Stress Resultant","$N_{ztheta}$: Shear Membrane Stress Resultant","$N_z$: Meridional Membrane Stress Resultant"];

N_zth_fig = plt.figure(1);
ax1 = N_zth_fig.add_subplot(projection="3d");
ax1.set_aspect("auto");
ax1.set_title(N_components[1]);
ax1.set(xticklabels=[],yticklabels=[]);
ax1.set_zlabel("z [mm]");
N_zth_max = 2*abs(T/(2*math.pi*(3000**2)) + Q/(math.pi*3000));

N_z_fig = plt.figure(2);
ax2 = N_z_fig.add_subplot(projection="3d");
ax2.set_aspect("auto");
ax2.set_title(N_components[2]);
ax2.set(xticklabels=[],yticklabels=[]);
ax2.set_zlabel("z [mm]");
N_z_max = 2*abs(M/(math.pi*(3000**2)) + P/(math.pi*2*3000));

# decision = input("Plot figures? Y/N ");

tol = 1e-3;

check = {ID: np.zeros((3,),dtype=bool) for i,ID in enumerate(strakeList)};
utilisation = {ID: np.empty((3,)) for i,ID in enumerate(strakeList)};
error = {ID:{"N_zth baseline":[],"N_zth variation":[],"N_z baseline":[],"N_z variation":[]} for k,ID in enumerate(strakeList)};
sigma_xRd = [];
tau_xthRd = [];
sigma_xEd = [];
tau_xthEd = [];

## loop over all strakes
for i,strakeID in enumerate(strakeList):
    print("Index "+str(i)+", ID= "+str(strakeID));

    # gives list of [H,R,t] values for given strake, in mm
    h, R, t = geometry.findStrakeGeometry(geoms_filename, strakeID);
    # global strake position
    z0 = geometry.findStrakePositionGlobal(geoms_filename, strakeID);

    ## LOCAL cylindrical coordinates r,th,z (radial, circumferential, meridional)
    Theta, Z = np.meshgrid(np.linspace(0,sec_angle,th_n), np.linspace(z0,z0+h,100),indexing="ij");
    X, Y = R*np.cos(Theta), R*np.sin(Theta); # for plotting only

    ## find global internal forces
    ### reaction forces R(th,z*=0) [N/mm]
    R_z = -(P/(sec_angle*R) + p_z*t*(H-z0))*np.ones(np.shape(Theta)); # vertical reaction force
    # R_x = - (Q - p_r*(H-z0)*np.cos(theta)); # p_r cancels over 2pi integral
    R_th = - (Q*np.sin(Theta)/(sec_angle*R) + p_th*h);
    R_r = - (Q*-np.cos(Theta)/(sec_angle*R) + p_r*h) ;

    R_T = -(T + p_th*sec_angle*R*R*h)*np.ones(np.shape(Theta)); # reaction torque [Nmm]
    R_M = -(M + Q*h)*np.ones(np.shape(Theta)); # reaction moment [Nmm]

    #### Nz0 i.e. N(z*=0) = -ve reaction forces
    Nz0 = np.array([-R_r,-R_th,-R_z,-R_M,-R_T]); # same order as "titles"

    ### internal forces at z*=h
    diagrams = {k : np.zeros(np.shape(Theta)) for k in titles};

    for k,key in enumerate(titles):
        diagrams[key] = Nz0[k,:,:]; # set z*=0 values correctly
        # also initialises other values to be edited below

    #### set z*=(0,h] values correctly
    diagrams[titles[0]] += -p_r*(Z-z0); # radial shear [N/mm]
    diagrams[titles[1]] += -p_th*(Z-z0); # circumferential shear [N/mm]
    diagrams[titles[2]] += -p_z*t*(h-Z-z0); # axial force (N/mm)
    diagrams[titles[3]] += -Q*(Z-z0); # bending moment (Nmm)
    diagrams[titles[4]] += -p_th*sec_angle*(R**2)*(Z-z0); # torque (Nmm)

    ## superimpose Ns
    ### from applied loading, N[0] = N_th, N[1] = N_zth, N[2] = N_z
    N = np.zeros(np.shape(Nz0))[0:3,:,:]; # initialise

    # adding fundamental results at all Z=[z0,h]
    N += fsr.p_rStresses(p_r, R, Z);
    N += fsr.p_thStresses(p_th,H-Z);
    N += fsr.p_zStresses(p_z + p_zSW*t, Z,H);
    N += fsr.PStresses(P, R, Z);
    N += fsr.TStresses(T, R, Theta);
    N += fsr.MStresses(M,R,Theta);
    N += fsr.QStresses(Q,R,Theta,H-Z);

    ## error between Ns and applied loads
    ### N_z
    error[strakeID]["N_z baseline"].append(abs(P - ((np.amax(N[2,:,-1])+np.amin(N[2,:,-1]))/2)*2*math.pi*R)/abs(P)); # base (not th dependent)
    error[strakeID]["N_z variation"].append(abs((abs(M) - abs(np.amax(N[2,:,-1])-np.amin(N[2,:,-1]))*0.5*math.pi*R*R)/M));
    ### N_zth
    if T != 0:
        error[strakeID]["N_zth baseline"].append(abs(T - ((np.amax(N[1,:,:])+np.amin(N[1,:,:]))/2)*sec_angle*R*R)/abs(T));
    error[strakeID]["N_zth variation"].append(abs(Q - (np.amax(N[1,:,:])-np.amin(N[1,:,:]))*0.5*math.pi*R)/abs(Q));

    tolcheck = 0;
    for k in error[strakeID].keys():
        tolcheck += sum(error[strakeID][k]);
    if tolcheck <= tol:
        print("Global Equilibrium satisfied.");
    else:
        print("Global equilibrium not satisfied for strake "+str(strakeID));
    # plotting

    # plot global loading "diagrams"
    for k,key in enumerate(titles):
        fig = plt.figure(num=k+3);
        ax = fig.add_subplot(projection="3d");
        ax.plot_surface(Theta,Z,diagrams[key]);

        # ax.set_title(key+", strake "+str(strakeID));

    # plot Ns
    # normalized colour distributions
    N_zth_cols = cm.jet((N[1,:,:]-np.amin(N[1,:,:]))/N_zth_max);
    N_z_cols = cm.jet((N[2,:,:]-np.amin(N[2,:,:]))/(2*N_z_max));

    ax1.plot_surface(X,Y,Z,facecolors=N_zth_cols,rstride=1, cstride=10);
    ax2.plot_surface(X,Y,Z,facecolors=N_z_cols,rstride=1, cstride=10);

    # buckling check LS3
    # axial
    sigma_xRd.append(LS3.findAxialBucklingStress(E, f_yk, Q_x["Class A"], lambda_x0, chi_xh, h, R, t, gamma_M1));
    sigma_xEd.append(np.amin(N[2,:,:])/t);
    check[strakeID][2], utilisation[strakeID][2] = LS3.checkIndividualStresses(sigma_xEd[i], sigma_xRd[i]);
    print("Axial buckling OK? "+str(check[strakeID][2]));
    print("Utilisation = "+str(utilisation[strakeID][2])+" %");

    # shear
    tau_xthRd.append(LS3.findShearBucklingStress(E, f_yk, Q_tau["Class A"], lambda_tau0, chi_tauh, beta_tau, eta_tau, h, R, t, gamma_M1));
    tau_xthEd.append(np.amin(N[1,:,:])/t);
    check[strakeID][1], utilisation[strakeID][1] = LS3.checkIndividualStresses(tau_xthEd[i], tau_xthRd[i]);
    print("Shear buckling OK? "+str(check[strakeID][1]));
    print("Utilisation = "+str(utilisation[strakeID][1])+" %\n");

    # stress interactions


# plot utilisation
fig = plt.figure();
ax = fig.add_subplot()
x = np.empty((2,len(strakeList)));
y = np.empty((2,len(strakeList)));
for i,ID in enumerate(strakeList):
    x[:,i] = [ID, utilisation[ID][1]];
    y[:,i] = [ID, utilisation[ID][2]];
ax.plot(x[0,:],x[1,:]);
ax.plot(y[0,:],y[1,:]);
ax.legend(["Shear","Axial"]);
ax.set_title("Utilisation [%]");
ax.set_xlabel("Strake ID");
# ax.set(xticklabels=x[0,:]);
plt.show();
