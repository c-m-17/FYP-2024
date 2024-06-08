"""
This script defines the strake class.

@author: Cerys Morley
"""
import conicalgeometrytocylindrical as geometry
import math

class strake:
    """
    Individual hollow cylindrical section of tower with constant thickness, height and radius.

    h: Height
    r: Average radius
    r_top: Top radius
    r_bot: Bottom radius
    t: Thickness
    beta: Apex half-angle
    z0: Global z-coordinate of the strake base (where z=0 is tower base).
    C_x: Critical buckling factor under axial compression
    C_tau: Critical buckling factor under shear.
    """
    def __init__(self,strakeID: str | int, fileName: str) -> None:
        """
        strakeID: unique ID for the strake
        fileName: relative path of file
        """
        self.h, d_1, d_2, self.t = geometry.findStrakeGeometry(fileName,strakeID)
        self.z0 : float = geometry.findStrakePositionGlobal(fileName,strakeID)

        self.r_top : float = d_1/2
        self.r_bot : float = d_2/2

        self.update()

    def update(self) -> None:
        """
        Updates strake's geometric properties when a new radius has been defined.
        """
        self.r : float = (self.r_top + self.r_bot)/2

    def calcC_x(self, omega: float) -> float:
        """
        Calculates critical buckling factor under axial compression, C_x.

        omega: Relative length of cylinder.

        Reference: EN 1993-1-6 Annex D.
        returns: self.C_x
        """
        if omega < 1.7: # EN 1993-1-6 D.3
            # short length
            C_x : float = 1.36 - (1.83/omega) + (2.07/(omega**2)) # EN 1993-1-6 D.8
        elif omega < 1.43*self.r/self.t: # EN 1993-1-6 D.4
            # medium length
            C_x : float = 1.0 # EN 1993-1-6 D.7  
        else:
            # long length
            C_x : float = 1.0
        return C_x

    def calcC_tau(self, omega: float) -> float:
        """
        Calculates critical buckling factor under torsion, C_tau.
        Assumes BC1r or BC2r conditions at the strake edges.

        omega: Relative length of cylinder.

        Reference: EN 1993-1-6 Annex D.
        """
        if omega < 10: # EN 1993-1-6 D.37
            # short cylinder
            # assuming BC1r or BC2r:
            alpha_taus : float = 120 - 130/(1+ 0.015*self.r/self.t) # EN 1993-1-6 D.43
            # assuming BC1f or BC2f:
            # alpha_taus : float = 70 - 75/(1 + 0.015*((r/t)**1.1)) # EN 1993-1-6 D.44

            b : float = 3 - 5/(1 + 0.4*((self.r/self.t)**0.6)) # EN 1993-1-6 D.45

            C_tau : float = math.sqrt(1 + alpha_taus/(omega**b)) # EN 1993-1-6 D.42

        elif omega < 8.7*self.r/self.t: # EN 1993-1-6 D.38
            # medium-length cylinder
            C_tau : float = 1.0 # EN 1993-1-6 D.41

        else: # EN 1993-1-6 D.39
            # long cylinder
            C_tau : float = math.sqrt(omega*self.t/self.r)/3 # EN 1993-1-6 D.46
        return C_tau

if __name__ == "__main__":  
    geoms_filename : str = "Sadowskietal2023-benchmarkgeometries.csv"

    strakeList, H = geometry.listStrakeIDs(geoms_filename)

    strakes : dict[str | int,strake] = {}
    for ID, strakeID in enumerate(strakeList):
        strakes[strakeID] = strake(strakeID,geoms_filename)

    print(strakes)