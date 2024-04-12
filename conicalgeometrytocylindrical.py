# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:16:51 2024

@author: Cerys Morley
"""
import math
import pandas as pd

def averageRadius(d_1,d_2):
    R = 0.25*(d_1+d_2)
    return R

# returns index & imported file of given strake
def findStrakeIndex(filename,strakeID):
    geoms = pd.read_csv(filename)
    I = geoms.loc[geoms["ID"] == strakeID].index.tolist()[0]
    return geoms, I

# def inclinedMeridionalLength(h,beta):
#     return h/math.cos(beta)

# returns list of strakeIDs from import file
def listStrakeIDs(filename):
    df = pd.read_csv(filename)
    StrakeIDlist = df["ID"]

    H = sum(df["h (mm)"]) # total tower height
    V = 0
    for i, ID in enumerate(StrakeIDlist):
        h, R, t = findStrakeGeometry(filename, ID)
        V += h*t*2*math.pi*R # tower volume

    return StrakeIDlist, H, V

# returns the *cylindrical* geometry of a given strake
def findStrakeGeometry(filename,strakeID):
    geoms, I = findStrakeIndex(filename, strakeID)

    h = geoms['h (mm)'].iloc[I]
    beta = geoms['beta (rad)'].iloc[I]
    d_1 = geoms['d1 (top) (mm)'].iloc[I]
    d_2 = geoms['d2 (bottom) (mm)'].iloc[I]
    t = geoms['t (mm)'].iloc[I]

    R = averageRadius(d_1, d_2)
    return [h,R,t]

# returns the global z-coordinate of the given strake's bottom boundary
def findStrakePositionGlobal(filename,strakeID):
    geoms, I = findStrakeIndex(filename, strakeID)

    H = sum(geoms["h (mm)"]) # total tower height
    z0 = H - sum(geoms["h (mm)"].iloc[range(0,I+1,1)]) # z*=0 in global z-coord

    return z0
