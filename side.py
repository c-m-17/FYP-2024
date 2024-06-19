"""
side.py

Uses main.py, LS3bucklingcheck.py and strake.py functions to calculate without optimisation.

Cerys Morley
2024
"""

# package imports
import numpy as np

# module imports
from main import *
from strake import strake
import LS3bucklingcheck as LS3

# set up input
nr = 3
nt = 3
radii = np.linspace(500,5500,nr)
thicknesses = np.linspace(2,50,nt)
geoms_filename : str = "restest-geometries.csv"
loads_filename : str = "loading-01.csv"

# initialise output
sigma_xRd = np.zeros((nr,nt))
tau_xthRd = sigma_xRd.copy()
V = sigma_xRd.copy()

# copied from main.py
# constants & parameters
f_yk : float = 355e6 # characteristic steel strength [N/m2]
gamma_M1 : float = 1.1 # material partial factor of safety (EN 1993-1-6 Table 4.2)
E : float = 210e9 # Young's Modulus [N/m2]
rho : float = 7850 # material density [kg/m3]
fabClass : str = "A" # fabrication quality class (one of ["A", "B", "C"])
beta_max : float = math.pi/8 # maximum apex half-angle [rad]

# import loading
loadNames : list[str] = ["p_r", "p_th", "p_z", "P", "Q", "T", "M"]
loads : dict[str,float] = {load : importLoadMagnitude(load,loads_filename) for load in loadNames}

for ir,r in enumerate(radii):
    for it,t in enumerate(thicknesses):
        # overwrite input geometry file with current r & t values
        with open(f"inputs/{geoms_filename}","w") as fil:
            fil.write(f"ID,h (mm),d1 (top) (mm),d2 (bottom) (mm),t (mm),\n")
            fil.write(f"0,2000,")
            fil.write(f"{r*2},{r*2},")
            fil.write(f"{t},\n")

        # import geometry
        strakeList = listStrakeIDs(geoms_filename)
        # create strake objects
        strakes : dict[str | int, strake] = {strakeID : strake(strakeID,geoms_filename) for i, strakeID in enumerate(strakeList)}

        # loop over strakes to find design resistance
        for s in strakes.values():
            sigma_xRd[ir][it] = LS3.findAxialBucklingStress(E,f_yk,fabClass,gamma_M1,s)[0]
            tau_xthRd[ir][it] = LS3.findShearBucklingStress(E,f_yk,fabClass,gamma_M1,s)[0]
            V[ir,it] = volumeFunction([s.t, s.r_bot],s)

            design = findMembraneStresses(s,loads,rho)
            print(design)


with open(f"sigma_xRd(rt).csv","w") as fil:
    for r in range(nr):
        for t in range(nt):
            fil.write(f"{sigma_xRd[r][t]},")
        fil.write(f"\n")

with open(f"tau_xthRd(rt).csv","w") as fil:
    for r in range(nr):
        for t in range(nt):
            fil.write(f"{tau_xthRd[r][t]},")
        fil.write(f"\n")

with open(f"volume(rt).csv","w") as fil:
    for r in range(nr):
        for t in range(nt):
            fil.write(f"{V[r][t]},")
        fil.write(f"\n")
