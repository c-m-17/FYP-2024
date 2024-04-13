"""
This script defines the strake class.

@author: Cerys Morley
Created on 19 Feb 2024
"""
import conicalgeometrytocylindrical as geometry

class strake:
    """Individual section of tower with constant thickness, height and radius."""
    def __init__(self,strakeID: int, geoms_filename: str):
        # height, radius, thickness [mm]
        self.h, self.r, self.t = geometry.findStrakeGeometry(geoms_filename,strakeID)
        # base z-coord in global system [mm]
        self.z0 = geometry.findStrakePositionGlobal(geoms_filename,strakeID)

# if file is run directly
if __name__ == "__main__":  
    geoms_filename = "Sadowskietal2023-benchmarkgeometries.csv"

    # strake list import
    strakeList, H, V = geometry.listStrakeIDs(geoms_filename)

    strakes = {}
    for ID, strakeID in enumerate(strakeList):
        # create strake in dictionary
        strakes[strakeID] = strake(strakeID)

    print(strakes)