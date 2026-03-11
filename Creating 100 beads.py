# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:14:20 2026

@author: charl
"""

import numpy as np
import matplotlib.pyplot as plt
import random

number_of_beads = 100

# Central arctic
x0=175
y0=155


def oil_beads(x0, y0, number_of_beads=100):
    #Set side length
    n = int(np.sqrt(number_of_beads))
    #Create 10 evenly spaced data points that are centred on 0 and offset from the edges by 0.5
    spacing = np.linspace(-0.5, 0.5, n, endpoint=False) + 0.5/n
    
    #Create a 10x10 square grid shape
    gx, gy = np.meshgrid(spacing, spacing)  
    
    #shape the grid about the coordinate points
    bead_x = x0 + gx.ravel()  # shape (100,)
    bead_y = y0 + gy.ravel()  # shape (100,)
    
    return bead_x, bead_y


bead_x, bead_y = oil_beads(x0=x0, y0=y0)
# print("X Values:")
# print(bead_x)
# print("-------------------------------")
# print("Y Values:")
# print(bead_y)
# print("-------------------------------")

plt.scatter(bead_x,bead_y,color='black',s=25)
plt.scatter(x0,y0,marker='*',s=200)
#plt.xlim(x0-5, x0+5)
#plt.ylim(y0-5,y0+5)
plt.title("Step number:")
plt.show()



