"""
This script defines the strake class.

@author: Cerys Morley
"""
import math
import pandas as pd

class strake:
    """
    Individual hollow cylindrical section of tower with constant thickness, height and radius.

    h: Height
    r: Average radius
    r_top: Top radius
    r_bot: Bottom radius
    t: Thickness
    L: Meridional length
    z0: Global z-coordinate of the strake base (where z=0 is tower base).
    C_x: Critical buckling factor under axial compression
    C_tau: Critical buckling factor under shear.
    """
    def __init__(self, strakeID: str | int, fileName: str) -> None:
        """
        strakeID: unique ID for the strake
        fileName: relative path of file
        """
        self.h, d_1, d_2, self.t, self.z0 = self.findStrakeGeometry(fileName,strakeID)

        self.update(r_top=d_1/2, r_bot=d_2/2)

    def update(self,*, r_top=None, r_bot=None) -> None:
        """
        Updates strake's geometric properties to edit the radii values.

        r_top: Radius of top of strake (float | int)
        r_bot: Radius of bottom of strake (float | int)
        """
        if r_top is not None:
            self.r_top = r_top

        if r_bot is not None:
            self.r_bot = r_bot
        
        self.r : float = (self.r_top + self.r_bot)/2

    def findStrakeGeometry(self, filename: str, strakeID: str | int) -> tuple[float]:
        """
        Approximates the geometry of a conical strake as cylindrical.

        filename: File name (in "inputs" folder)
        strakeID: unique ID of strake

        returns: Height, Top diameter, Bottom diameter, Thickness (all in metres), z0 position
        """
        # dataframe & corresponding strake index
        geoms : pd.DataFrame = pd.read_csv(f"inputs\{filename}")
        I : str | int = geoms.loc[geoms["ID"] == strakeID].index.tolist()[0]

        h : float = geoms['h (mm)'].iloc[I]/1000
        d_1 : float = geoms['d1 (top) (mm)'].iloc[I]/1000
        d_2 : float = geoms['d2 (bottom) (mm)'].iloc[I]/1000
        t : float = geoms['t (mm)'].iloc[I]/1000

        H : float = sum(geoms["h (mm)"]) # total tower height
        z0 : float = (H - sum(geoms["h (mm)"].iloc[range(0,I+1,1)]))/1e3 # z*=0 in global z-coord

        return h, d_1, d_2, t, z0

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
