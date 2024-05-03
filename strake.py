"""
This script defines the strake class.

@author: Cerys Morley
"""
import conicalgeometrytocylindrical as geometry
import math

class strake:
    """
    Individual hollow cylindrical section of tower with constant thickness, height and radius.

    h: height
    r: radius
    t: thickness
    z0: Global z-coordinate of the strake base (where z=0 is tower base).
    omega: relative length of strake
    C_x: critical buckling factor under axial compression
    C_tau: critical buckling factor under shear.
    """
    def __init__(self,strakeID: str | int, fileName: str) -> None:
        """
        strakeID: unique ID for the strake
        fileName: relative path of file
        """
        self.h, self.r, self.t = geometry.findStrakeGeometry(fileName,strakeID)
        self.z0 = geometry.findStrakePositionGlobal(fileName,strakeID)

        self.omega = self.h/math.sqrt(self.r*self.t) # EN 1993-1-6 D.1

        self.C_x = self.calcC_x()
        self.C_tau = self.calcC_tau()

    def calcC_x(self) -> float:
        """
        Calculates critical buckling factor under axial compression, C_x.

        Reference: EN 1993-1-6 Annex D.
        """
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

    def calcC_tau(self) -> float:
        """
        Calculates critical buckling factor under shear, C_tau.
        Assumes BC1r or BC2r conditions at the strake edges.

        Reference: EN 1993-1-6 Annex D.
        """
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

if __name__ == "__main__":  
    geoms_filename : str = "Sadowskietal2023-benchmarkgeometries.csv"

    strakeList, H, V = geometry.listStrakeIDs(geoms_filename)

    strakes : dict[str | int,strake] = {}
    for ID, strakeID in enumerate(strakeList):
        strakes[strakeID] = strake(strakeID,geoms_filename)

    print(strakes)