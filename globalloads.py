# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:20:40 2024

@author: Cerys Morley
"""
# NOTE: this function can only handle constant load values!!!

import pandas as pd

def importLoadMagnitude(loadtype,filename):
    
    loads = pd.read_csv(filename);
    # loadtype is string
    A = loads["Magnitude"].loc[loads["Load"] == loadtype];
    
    return A.unique()[0]
