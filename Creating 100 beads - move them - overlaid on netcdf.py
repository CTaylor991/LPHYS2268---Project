# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:14:20 2026

@author: charl
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import xarray as xr



number_of_beads = 100

sea_ice_min = 0.15

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
    bead_x = x0 + gx.ravel()  
    bead_y = y0 + gy.ravel()  
    
    return bead_x, bead_y

#-------Loading in NetCDF fil and plotting for background--------------------------------------
ds = xr.open_dataset("siconc_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")
siconc = ds['siconc']   # shape (366, 361, 361)
#lat = ds['lat'].values  # shape (361, 361)
#lon = ds['lon'].values  # shape (361, 361)

#----------Build custom colourmap-----------------------------
cmap = plt.cm.Blues_r
#land colours defined as nan
cmap.set_bad(color='darkgreen')   
#Set colour for when threshold is broken (open ocean, no sea ice)
cmap.set_under(color='midnightblue')




#-------Bead creation and manipialation-----------------------------
bead_x, bead_y = oil_beads(x0=x0, y0=y0)
# print("X Values:")
    # print(bead_x)
    # print("-------------------------------")
    # print("Y Values:")
    # print(bead_y)
    # print("-------------------------------")
fig, ax = plt.subplots()
siconc.isel(time_counter=0).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
ax.invert_yaxis()

ax.scatter(bead_x,bead_y,color='black',s=25)
ax.scatter(x0,y0,marker='*',color='sienna',s=100)
#ax.set_xlim(x0-1, x0+30)
#ax.set_ylim(y0-1, y0+30)
ax.set_title("Starting condition")
plt.show()

#---------------------------------------------------------------
steps = 360

#initialise arrays to store movement rows are time steps and xolumns are beads
movement_x = np.zeros((steps + 1, number_of_beads))
movement_y = np.zeros((steps + 1, number_of_beads))

#Store the initial position
movement_x[0] = bead_x
movement_y[0] = bead_y

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
        # fin_x = bead_x[i]
        # fin_y = bead_y[i]
        if i == 0:
            if (j+1) % 25 ==0:
                print("")
                print("Max change in x for bead",i+1," for step ", j+1 , "=", 0.65*(init_x/x0))
                print("Max Change in y for bead",i+1," for step ", j+1 , "=", 0.65*(init_y/y0))
    #Save the new position to the trajectory arrays. +1 because it is the next step in the array and to preserve x0 and y0
    movement_x[j+1] = bead_x
    movement_y[j+1] = bead_y
    
    if (j+1) % 36 == 0:
        fig, ax = plt.subplots()
        
        # Background sea ice for this timestep
        siconc.isel(time_counter=j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
        ax.invert_yaxis()
        
        #Plot the new position of the beads and the centre point
        ax.scatter(x0,y0,marker='*',s=100,color='sienna',zorder=3)    
        ax.scatter(bead_x,bead_y,color='black',s=25,zorder=2)
    
        #Plotting of moevement of oil beads
        for i in range(number_of_beads):
            #Plot the trajectories for each bead i being bead and j being steps +2 to include all steps up to the last
            ax.plot(movement_x[:j+2, i], movement_y[:j+2, i], alpha=0.8, linewidth=0.5,zorder=1)
            
        #ax.set_xlim(x0-1, x0+30)
        #ax.set_ylim(y0-1, y0+30)
        ax.set_title(f"Step number:{j+1}")
        plt.show()

print(siconc.isel(time_counter=0).values.max())
print(siconc.isel(time_counter=0).values.min())
print(np.unique(siconc.isel(time_counter=0).values[:10,:10]))  # sample corner values


