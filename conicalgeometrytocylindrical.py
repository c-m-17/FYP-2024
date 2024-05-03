# -*- coding: utf-8 -*-
"""
Functions to extract the geometry of any number of truncated conical strakes in a tower,
and approximate their cylindrical equivalents.
Also extracts loads from files.

@author: Cerys Morley
"""
import math
import pandas as pd

def averageRadius(d_top: float, d_bottom: float) -> float:
    """
    Calculates average strake radius, assuming linear variation over the strake height.

    d_top: top diameter of strake
    d_bottom: bottom diameter of strake
    """
    R = 0.25*(d_top+d_bottom)

    return R

def findStrakeIndex(filename: str, strakeID: str | int) -> tuple[pd.DataFrame,int]:
    """
    Finds index of strake within a dataframe, and returns it along with the file
    itself as a DataFrame.

    filename: relative file path
    strakeID: ID from "ID" column in dataframe for given strake

    returns: strake geometry DataFrame, index of given strake in DataFrame
    """
    geoms : pd.DataFrame = pd.read_csv(filename)
    I : int = geoms.loc[geoms["ID"] == strakeID].index.tolist()[0]

    return geoms, I

def listStrakeIDs(filename: str) -> tuple[pd.Series,float,float]:
    """
    Extracts strake IDs from file. Calculates tower height and volume.

    filename: relative file path

    returns: strake IDs, tower height, tower volume
    """
    df : pd.DataFrame = pd.read_csv(filename)
    StrakeIDlist : pd.Series = df["ID"]

    H : float = sum(df["h (mm)"])/1e3 # total tower height
    V : float = 0

    for i, ID in enumerate(StrakeIDlist):
        h, R, t = findStrakeGeometry(filename, ID)
        V += h*t*2*math.pi*R # tower volume

    return StrakeIDlist, H, V

def findStrakeGeometry(filename: str, strakeID: str | int) -> tuple[float,float,float]:
    """
    Approximates the geometry of a conical strake as cylindrical.

    filename: relative file path
    strakeID: unique ID of strake

    returns: Height, Radius, Thickness
    """
    geoms, I = findStrakeIndex(filename, strakeID)

    h : float = geoms['h (mm)'].iloc[I]/1000
    beta : float = geoms['beta (rad)'].iloc[I]
    d_1 : float = geoms['d1 (top) (mm)'].iloc[I]
    d_2 : float = geoms['d2 (bottom) (mm)'].iloc[I]
    t : float = geoms['t (mm)'].iloc[I]/1000

    R : float = averageRadius(d_1, d_2)/1000
    return h, R, t

# returns the global z-coordinate of the given strake's bottom boundary
def findStrakePositionGlobal(filename: str, strakeID: str | int) -> float:
    """
    Finds global z-coordinate of the bottom of a strake, where z=0 is at base of tower.
    
    filename: relative file path
    strakeID: unique ID of strake

    returns: z-coordinate (vertical) of strake base
    """
    geoms, I = findStrakeIndex(filename, strakeID)

    H : float = sum(geoms["h (mm)"]) # total tower height
    z0 : float = (H - sum(geoms["h (mm)"].iloc[range(0,I+1,1)]))/1e3 # z*=0 in global z-coord

    return z0

def importLoadMagnitude(loadType: str, filename: str) -> float:
    """
    Extracts the magnitude of a specific load type from a file.
    File must have "Load" & "Magnitude" columns.

    loadType: One of [p_r, p_th, p_z, P, T, Q, M]
    filename: Relative file path

    returns: Global constant value of that load.
    """
    
    loads : pd.DataFrame = pd.read_csv(filename)
    A = loads["Magnitude"].loc[loads["Load"] == loadType]
    
    return A.unique()[0]