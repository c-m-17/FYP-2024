# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:16:51 2024

@author: Cerys Morley
"""
import math
import pandas as pd

def averageRadius(d_1,d_2,beta):
    return 0.25*(d_1+d_2)/math.cos(beta)

def inclinedMeridionalLength(h,beta):
    return h/math.cos(beta)

# returns the *cylindrical* geometry of a given strake
def findStrakeGeometry(filename,strakeID):
    geoms = pd.read_csv(filename);

    h = geoms['h (mm)'].loc[geoms["ID"] == strakeID];
    beta = geoms['beta (rad)'].loc[geoms["ID"] == strakeID];
    d_1 = geoms['d1 (top) (mm)'].loc[geoms["ID"] == strakeID];
    d_2 = geoms['d2 (bottom) (mm)'].loc[geoms["ID"] == strakeID];
    t = geoms['t (mm)'].loc[geoms["ID"] == strakeID];
    
    l = inclinedMeridionalLength(h, beta);
    R = averageRadius(d_1, d_2, beta);
    return [l.unique(),R.unique(),t.unique()]
