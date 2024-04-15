# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 16:57:35 2024

@author: Cerys Morley

This file creates a generic plot with formatting variables as inputs.
"""

import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()

def plotter(figno,ncols,nrows,title,x,xlabel,y,ylabel):
    fig, ax = plt.subplots(num = figno,layout="constrained")
    
    ax.plot(x,y);
    
    ax.set_title(title);
    ax.set_xlabel(xlabel);
    ax.set_ylabel(ylabel);
    ax.legend();
    
    return fig

# NOTE currently only works with nrows = 1 and ncols = 1.