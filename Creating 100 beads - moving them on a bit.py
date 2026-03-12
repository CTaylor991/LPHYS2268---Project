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


plt.scatter(bead_x,bead_y,color='black',s=25)
plt.scatter(x0,y0,marker='*',s=200)
plt.xlim(x0-1, x0+30)
plt.ylim(y0-1, y0+30)
plt.title("Step number:")
plt.show()

#---------------------------------------------------------------
#Accessing individual beads and their data and then changing the position
steps = 100

#initialise arrays to store trajectories rows are time steps and xolumns are beads
traj_x = np.zeros((steps + 1, number_of_beads))
traj_y = np.zeros((steps + 1, number_of_beads))

#Store the initial position --CHARLIE TO CHECK THIS IS THE GRIDDED POINTS
traj_x[0] = bead_x
traj_y[0] = bead_y

for j in range(steps):
    for i in range(number_of_beads):
        #print("Bead",i+1,"location:",bead_x[i],bead_y[i])
        
        #add a random value to the x and y position and save as the new position
        #bead_x[i] = bead_x[i] + random.uniform(-0.25, 0.5)
        #bead_y[i] = bead_y[i] + random.uniform(-0.25, 0.5)
        init_x = bead_x[i]
        init_y = bead_y[i]
        bead_x[i] = bead_x[i] + random.uniform(-0.25*(bead_x[i]/x0), 0.65*(bead_x[i]/x0))
        bead_y[i] = bead_y[i] + random.uniform(-0.25*(bead_y[i]/y0), 0.65*(bead_y[i]/y0))
        fin_x = bead_x[i]
        fin_y = bead_y[i]
        if i == 0:
            if (j+1) % 25 ==0:
                print("")
                print("Max change in x for bead",i+1," for step ", j+1 , "=", 0.65*(init_x/x0))
                print("Max Change in y for bead",i+1," for step ", j+1 , "=", 0.65*(init_y/y0))
    #Save the new position to the trajectory arrays. +1 because it is the next step in the array and to preserve x0 and y0
    traj_x[j+1] = bead_x
    traj_y[j+1] = bead_y
    
    if (j+1) % 10 == 0:
        #Plot the new position of the beads and the centre point
        plt.scatter(x0,y0,marker='*',s=200,zorder=3)    
        plt.scatter(bead_x,bead_y,color='black',s=25,zorder=2)
    
        #Plotting of trajectories
        for i in range(number_of_beads):
            #Plot the trajectories for each bead i being bead and j being steps +2 to include all steps up to the last
            plt.plot(traj_x[:j+2, i], traj_y[:j+2, i], alpha=0.8, linewidth=0.5,zorder=1)
            
        plt.xlim(x0-1, x0+30)
        plt.ylim(y0-1, y0+30)
        plt.title(f"Step number:{j+1}")
        plt.show()



