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
loads_filename = "Sadowskietal2023-benchmarkloads-LC2.csv";

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
f_yk = 355; # characteristic steel strength [N/mm2]
gamma_M1 = 1.1; # EN 1993-1-6 Table 4.2
E = 210e3; # Young's Modulus [N/mm2]

Q_x = {"Class A":40, "Class B":25, "Class C": 16}; # EN 1993-1-6 Table D.1
chi_xh = 1.10; # EN 1993-1-6 D.19
lambda_x0 = 0.10; # EN 1993-1-6 D.10

# find number of strakes
strakeList, H = geometry.listStrakeIDs(geoms_filename);

titles = ["Shear (radial) Force","Shear (theta) Force","Axial Force Diagram","Bending Moment Diagram","Global Torque Diagram"];
N_components = ["$N_{\theta}$: Circumferential Membrane Stress Resultant","$N_{z\theta}$: Shear Membrane Stress Resultant","$N_z$: Meridional Membrane Stress Resultant"];

tower_fig = plt.figure(1);
decision = input("Plot figures? Y/N ");

## loop over all strakes
for i,strakeID in enumerate(strakeList):
    print("Index "+str(i)+", ID= "+str(strakeID));

    # gives list of [H,R,t] values for given strake, in mm
    h, R, t = geometry.findStrakeGeometry(geoms_filename, strakeID);
    # global strake position
    z0 = geometry.findStrakePositionGlobal(geoms_filename, strakeID);

    ## LOCAL cylindrical coordinates r,th,z (radial, circumferential, meridional)
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
    diagrams[titles[2]] += -p_z*(Z-z0); # axial force (N/mm)
    diagrams[titles[3]] += -Q*(Z-z0); # bending moment (Nmm)
    diagrams[titles[4]] += -p_th*sec_angle*(R**2)*(Z-z0); # torque (Nmm)

    ## superimpose Ns
    ### from applied loading, N[0] = N_th, N[1] = N_zth, N[2] = N_z
    N = np.zeros(np.shape(Nz0))[0:3,:,:]; # initialise

    # adding fundamental results at all Z=[z0,h]
    N += fsr.p_rStresses(p_r, R, Z);
    N += fsr.p_thStresses(p_th, Z);
    N += fsr.p_zStresses(p_z, Z);
    N += fsr.PStresses(P, R, Z);
    N += fsr.TStresses(T, R, Theta);
    N += fsr.MStresses(M,R,Theta);
    N += fsr.QStresses(Q,R,Theta,h-Z);

    # output values at z*=H
    print("\nmax N_th(z*=h) = "+str(max(N[0,:,-1])));
    print("\nmax N_zth(z*=h) = "+str(max(N[1,:,-1])));
    print("\nmax N_z(z*=h) = "+str(max(N[2,:,-1])));

    ## error between Ns and diagrams
    tol = 1e-3;
    print("Global Equilibrium satisfied?\n");
    ### N_z
    absrelerr = abs(P - (np.mean(N[2,:,:])*2*math.pi*R))/abs(P); # base (not th dependent)
    print("Error for N_z baseline = " + str(absrelerr));
    print(absrelerr<=tol);
    absrelerr = abs(M - (np.amax(N[2,:,-1])-np.amin(N[2,:,-1]))*0.5*math.pi*R*R)/abs(M);
    print("Error for N_z variation = " + str(absrelerr));
    print(absrelerr<=tol);
    ### N_zth
    absrelerr = abs(T - (np.mean(N[1,:,:])*sec_angle*R*R))/abs(T);
    print("Error for N_zth baseline = " + str(absrelerr));
    print(absrelerr<=tol);
    absrelerr = abs(Q - (np.amax(N[1,:,:])-np.amin(N[1,:,:]))*0.5*math.pi*R)/abs(Q);
    print("Error for N_zth variation = " + str(absrelerr));
    print(absrelerr<=tol);

    # plotting

    ## plot strake geometry
    ax = tower_fig.add_subplot(projection="3d");
    ax.plot_surface(X,Y,Z);
    ax.set_aspect('auto');

    ax.set(xticklabels=[],
           yticklabels=[]);
    # ax.annotate("Height "+str(h)+" mm,\nRadius "+str(R)+" mm",(50,200),xycoords="figure points");
    plt.show();

    if decision == "Y":
        # plot global loading "diagrams"
        for k,key in enumerate(titles):
            fig = plt.figure();
            ax = fig.add_subplot(projection="3d");
            ax.plot_surface(Theta,Z,diagrams[key]);

            ax.set_title(key+", strake "+str(strakeID));
            plt.show();

        # plot Ns
        for i in range(3):
            fig = plt.figure();
            ax = fig.add_subplot(projection="3d")

            ax.plot_surface(Theta,Z,N[i,:,:],cmap=cm.coolwarm);
            ax.set_title(N_components[i]+", strake "+str(strakeID));

            ax.annotate("Max: "+str(np.amax(N[i,:,:]))+" N/mm\nMin: "+str(np.amin(N[i,:,:]))+" N/mm",(50,200),xycoords="figure points");
            plt.show();

        # fig = plt.figure();
        # ax = fig.add_subplot(projection="3d")

        # ax.plot_surface(Theta,Z,N[2,:,:],cmap=cm.coolwarm);
        # ax.set_title(N_components[2]+", Strake "+str(strakeID));

        # ax.annotate("Max: "+str(np.amax(N[2,:,:]))+" N/mm\nMin: "+str(np.amin(N[2,:,:]))+" N/mm",(50,200),xycoords="figure points");
        # plt.show();

    # buckling check LS3
    check, utilisation = LS3.checkAxialBuckling(E, f_yk, Q_x, lambda_x0, chi_xh, h, R, t);
    print("Axial buckling OK? "+str(check));
    print("Utilisation = "+str(utilisation)+" %\n");
