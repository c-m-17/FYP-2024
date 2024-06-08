# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 16:57:35 2024

@author: Cerys Morley

This file creates a generic plot with formatting variables as inputs.
"""

from matplotlib.axes import Axes
from matplotlib.figure import Figure
import seaborn as sns

sns.set_theme(style="whitegrid")

def plotTowerGeometry(fig : Figure, ax : Axes, strakes : dict, title : str) -> None:
    """
    Plots a single figure of the tower's geometry, and saves to a png file.

    fig: Fugure object
    ax: Axes object to plot to
    strakes: strakeID, strake object key-value pairs.
    title: Axes title

    returns: Figure of the rho-z plane.
    """
    # fig, (ax1, ax2) = plt.subplots(num = figno, ncols=2, sharey=True, layout="constrained")

    for ID, s in strakes.items():
        # strake joints
        ax.hlines(s.z0, -s.r_bot/2, s.r_bot/2,"blue")
        ax.hlines(s.z0 + s.h, -s.r_top/2, s.r_top/2,"blue")

        # strake sides
        ax.plot([-s.r_bot/2, -s.r_top/2], [s.z0, s.z0 + s.h], "b-",
                [s.r_bot/2, s.r_top/2], [s.z0, s.z0 + s.h], "b-")

    ax.set_title(f"{title} Support Tower Geometry")
    ax.set_xlabel(f"Strake radius [m]")

    sns.despine(fig,ax,trim=True,offset=10)
    fig.savefig(f"tower-geometry.png", format="png")