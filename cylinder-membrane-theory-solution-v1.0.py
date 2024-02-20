# -*- coding: utf-8 -*-
"""
This script finds the membrane stress resultants N acting on a cylindrical shell.
The geometry, applied loading and material properties are taken as inputs.
Only *membrane* equilibrium is considered.

author: Cerys Morley
date: 19 Feb 2024
"""
## import packages
import numpy as np

## input cylindrical shell geometry
# - radius "R" (mm)
# - height "H" (mm)
# - sector angle "T" (radians)
R = 500;
H = 1000;
T = 2*np.pi;
t = 1;

## input applied loading
# in cylindrical coordinates r,th,z (radial, circumferential, meridional)
# loads = f(th,z) *only*, no variation in r considered.
# arrays where th = row, z = column
th_n = 360; z_n = int(H/2)+1;
th = np.linspace(0,T,th_n);
z = np.linspace(0,H,z_n);
# coords[0] = z-coords, coords[1] = th-coords
coords = np.meshgrid(z,th);

# - pressures p_r, p_th, p_z (N/mm2 (= MPa))
p_r = (1 + 1*(np.divide(coords[0],H)))
p_th = 0*coords[0];
p_z = -0.5*p_r;

# - axial force P (N)
P = 0; # (+ve z-dir)

# - uniform torque T (Nmm (= kNm))
T = 0; # (+ve th direction)

# - global bending moment M (Nmm (=kNm))
M = 0; # (AC)

# - transverse shear force Q (N)
Q = 0; # (+ve x-dir)

## find global loading (superposition)
# reaction forces R(th)
R_z = -P/(2*R*np.pi) - p_z*H*R; # N
R_x = -Q; # N
R_M = M-Q*H; # Nmm (C)

# internal force distributions
N = -Q; # axial force (N)
V = P; # shear force (N)
m = -M-Q*(H-z); # bending moment (Nmm) 

## plot global loading
# - use matplotlib
# - function to add standard plot formatting
# - plot AFD on z-axis
# - plot SFD on z-axis
# - plot BMD on z-axis
# - plot torque on r-t plane?

## membrane equilibrium
# - radial equilibrium
N_th = p_r*R;

# - circumferential equilbrium
dp_rdth = np.gradient(p_r,axis=0);
N_zth = -np.trapz(dp_rdth+p_th,axis=1);
N_zth = np.transpose(N_zth*np.ones((z_n,th_n)));
f_1 = -N_zth; # N_zt(z=H) = 0 (BC1)
N_zth += f_1;

# - meridional equilibrium
dN_zthdth = np.gradient(N_zth,axis=0);
N_z = -np.trapz(dN_zthdth/R + p_z,axis=1)
f_2 = -N_z; # N_z(z=H) = 0 (BC2)
N_z += f_2;

# plot N distribution
# - plot on z-axis

# find critical N
# - max()

# identify failure mechanism
# - bursting or buckling?
# - location