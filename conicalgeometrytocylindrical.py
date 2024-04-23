# -*- coding: utf-8 -*-
"""
Extracts geometry of strakes in tower.

@author: Cerys Morley
"""
import math
import pandas as pd

def averageRadius(d_1: float, d_2: float):
    """Calculates average strake radius."""
    R = 0.25*(d_1+d_2)
    return R

# returns index & imported file of given strake
def findStrakeIndex(filename: str, strakeID: int):
    """Finds index of strake within dataframe."""
    geoms = pd.read_csv(filename)
    I = geoms.loc[geoms["ID"] == strakeID].index.tolist()[0]
    return geoms, I

# def inclinedMeridionalLength(h,beta):
#     return h/math.cos(beta)

# returns list of strakeIDs from import file
def listStrakeIDs(filename: str):
    "Creates pandas series of strake IDs"
    df = pd.read_csv(filename)
    StrakeIDlist = df["ID"]

    H = sum(df["h (mm)"])/1e3 # total tower height
    V = 0
    for i, ID in enumerate(StrakeIDlist):
        h, R, t = findStrakeGeometry(filename, ID)
        V += h*t*2*math.pi*R # tower volume

    return StrakeIDlist, H, V

# returns the *cylindrical* geometry of a given strake
def findStrakeGeometry(filename: str, strakeID: int):
    """Finds geometry of a strake."""
    geoms, I = findStrakeIndex(filename, strakeID)

    h = geoms['h (mm)'].iloc[I]/1000
    beta = geoms['beta (rad)'].iloc[I]
    d_1 = geoms['d1 (top) (mm)'].iloc[I]
    d_2 = geoms['d2 (bottom) (mm)'].iloc[I]
    t = geoms['t (mm)'].iloc[I]/1000

    R = averageRadius(d_1, d_2)/1000
    return h, R, t

# returns the global z-coordinate of the given strake's bottom boundary
def findStrakePositionGlobal(filename: str, strakeID: int):
    """Finds global z-coordinate of strake base"""
    geoms, I = findStrakeIndex(filename, strakeID)

    H = sum(geoms["h (mm)"]) # total tower height
    z0 = (H - sum(geoms["h (mm)"].iloc[range(0,I+1,1)]))/1e3 # z*=0 in global z-coord

    return z0
