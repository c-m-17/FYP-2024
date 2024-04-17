"""
This script defines the strake class.

@author: Cerys Morley
Created on 19 Feb 2024
"""
import conicalgeometrytocylindrical as geometry
import math

class strake:
    """Individual section of tower with constant thickness, height and radius."""
    def __init__(self,strakeID: int, geoms_filename: str):
        # height, radius, thickness [mm]
        self.h, self.r, self.t = geometry.findStrakeGeometry(geoms_filename,strakeID)
        # base z-coord in global system [mm]
        self.z0 = geometry.findStrakePositionGlobal(geoms_filename,strakeID)

        self.omega = self.h/math.sqrt(self.r*self.t)

        self.C_x = self.calcC_x()
        self.C_tau = self.calcC_tau()

        self.radiusConstrained = False

    def calcC_x(self):
        if self.omega < 1.7:# EN 1993-1-6 D.3
            # short length
            C_x = 1.36 - (1.83/self.omega) + (2.07/(self.omega**2)) # EN 1993-1-6 D.8
        elif self.omega < 1.43*self.r/self.t: # EN 1993-1-6 D.4
            # medium length
            C_x = 1 # EN 1993-1-6 D.7  
        else:
            # long length
            C_x = 1
        return C_x

    def calcC_tau(self):
        if self.omega < 10: # EN 1993-1-6 D.37
            # short cylinder
            # assuming BC1r or BC2r:
            alpha_taus = 120 - 130/(1+ 0.015*self.r/self.t) # EN 1993-1-6 D.43
            # assuming BC1f or BC2f:
            # alpha_taus = 70 - 75/(1 + 0.015*((r/t)**1.1)) # EN 1993-1-6 D.44

            b = 3 - 5/(1 + 0.4*((self.r/self.t)**0.6)) # EN 1993-1-6 D.45
            C_tau = math.sqrt(1 + alpha_taus/(self.omega**b)) # EN 1993-1-6 D.42
        elif self.omega < 8.7*self.r/self.t: # EN 1993-1-6 D.38
            # medium-length cylinder
            C_tau = 1.0 # EN 1993-1-6 D.41
        else: # EN 1993-1-6 D.39
            # long cylinder
            C_tau = math.sqrt(self.omega*self.t/self.r)/3 # EN 1993-1-6 D.46
        return C_tau

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