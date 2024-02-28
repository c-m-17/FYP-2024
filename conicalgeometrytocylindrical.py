# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:16:51 2024

@author: Cerys Morley
"""
import math
import pandas as pd

def averageRadius(d_1,d_2,beta):
    return 0.25*(d_1+d_2)/math.cos(beta)

# def inclinedMeridionalLength(h,beta):
#     return h/math.cos(beta)

# returns the *cylindrical* geometry of a given strake
def findStrakeGeometry(filename,strakeID):
    geoms = pd.read_csv(filename);
    I = geoms.loc[geoms["ID"] == strakeID].index.tolist()[0];
    
    h = geoms['h (mm)'].iloc[I];
    beta = geoms['beta (rad)'].iloc[I];
    d_1 = geoms['d1 (top) (mm)'].iloc[I];
    d_2 = geoms['d2 (bottom) (mm)'].iloc[I];
    t = geoms['t (mm)'].iloc[I];
    
    R = averageRadius(d_1, d_2, beta);
    return [h,R,t]

# returns the global z-coordinate of the given strake's bottom boundary
def findStrakePositionGlobal(filename,strakeID):
    geoms = pd.read_csv(filename);
    I = geoms.loc[geoms["ID"] == strakeID].index.tolist()[0];
    
    H = sum(geoms["h (mm)"]); # total tower height
    z0 = H - sum(geoms["h (mm)"].iloc[range(0,I+1,1)]); # z*=0 in global z-coord
    
    return [z0,H]
