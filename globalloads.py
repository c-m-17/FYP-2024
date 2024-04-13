"""
To import constant load values into a dictionary.

@author: Cerys Morley
"""
import pandas as pd

def importLoadMagnitude(loadtype: str, filename: str):
    
    loads = pd.read_csv(filename)
    # loadtype is string
    A = loads["Magnitude"].loc[loads["Load"] == loadtype]
    
    return A.unique()[0]

def allLoads(loads_filename: str, keys: list[str]):
    loads = {key : importLoadMagnitude(key,loads_filename) for k, key in enumerate(keys)}
    return loads
