# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:14:20 2026

@author: charl
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import xarray as xr

#Variables------------------------------------------
#Note 1st of November is j = 153
start_date = 153

number_of_beads = 100

steps = 365-start_date

sea_ice_min = 0.15

# Central arctic
x0=175
y0=155
# x0=200
# y0=200

#Create bead formation------------------------
#Set side length
n = int(np.sqrt(number_of_beads))
#Create 10 evenly spaced data points that are centred on 0 and offset from the edges by 0.5
spacing = np.linspace(-0.5, 0.5, n, endpoint=False) + 0.5/n
    
#Create a 10x10 square grid shape
gx, gy = np.meshgrid(spacing, spacing)  
    
#shape the grid about the coordinate points
bead_x = x0 + gx.ravel()  
bead_y = y0 + gy.ravel()  

#-------Loading in NetCDF files and plotting for background--------------------------------------
ds_i = xr.open_dataset("siconc_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")
ds_u = xr.open_dataset("sivelu_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")
ds_v = xr.open_dataset("sivelv_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")

siconc = ds_i['siconc']   

#----------Build custom colourmap-----------------------------
cmap = plt.cm.Blues_r
#land colours defined as nan
cmap.set_bad(color='darkgreen')   
#Set colour for when threshold is broken (open ocean, no sea ice)
cmap.set_under(color='midnightblue')


#Set up the background plot------------------------------------
fig, ax = plt.subplots()
siconc.isel(time_counter=start_date).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
ax.invert_yaxis()

ax.scatter(bead_x,bead_y,color='black',s=25)
ax.scatter(x0,y0,marker='*',color='sienna',s=100)
#ax.set_title("Starting condition")
plt.show()


#initialise arrays to store movement rows are time steps and xolumns are beads
movement_x = np.zeros((steps + 1, number_of_beads))
movement_y = np.zeros((steps + 1, number_of_beads))

#initialise arrays to store ice concentration rows are time steps and xolumns are ice con
iceconc_array = np.zeros((steps + 1, number_of_beads))

#Store the initial position
movement_x[0] = bead_x
movement_y[0] = bead_y

for i in range(number_of_beads):
    grid_x = int(bead_x[i])
    grid_y = int(bead_y[i])
    iceconc_array[0, i] = siconc.isel(time_counter=start_date, y=grid_y, x=grid_x).values


for j in range(steps):
    for i in range(number_of_beads):
        init_x = bead_x[i]
        init_y = bead_y[i]

                
        grid_x = int(init_x)
        grid_y = int(init_y)
        
        u_velocity = ds_u["sivelu"]
        u_velocity = u_velocity.isel(time_counter=start_date+j, y=grid_y, x=grid_x)
        
        v_velocity = ds_v["sivelv"]
        v_velocity = v_velocity.isel(time_counter=start_date+j, y=grid_y, x=grid_x)
        
        ice_conc = siconc.isel(time_counter=start_date+j, y=grid_y, x=grid_x)
        iceconc_array[j+1, i] = ice_conc.values
        if ice_conc < sea_ice_min:
            print("Bead ", i+1,"dropped out at step", j+1)
        init_x = bead_x[i]
        init_y = bead_y[i]
        #Move along (86,400 seconds in a day, 25km grid size)
        bead_x[i] = bead_x[i] + (u_velocity.values * 86400/25000)
        bead_y[i] = bead_y[i] + (v_velocity.values * 86400/25000)
    
               
        
        
    #Save the new position to the trajectory arrays. +1 because it is the next step in the array and to preserve x0 and y0
    movement_x[j+1] = bead_x
    movement_y[j+1] = bead_y
    
    
    
    if (j+1) % 20 == 0:

        fig, ax = plt.subplots()
        
        # Background sea ice for this timestep
        siconc.isel(time_counter=start_date+j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
        
        
        #Plot the new position of the beads and the centre point
        ax.scatter(x0,y0,marker='*',s=100,color='sienna',zorder=3)    
        ax.scatter(bead_x,bead_y,color='black',s=25,zorder=2)
        
    
        #Plotting of moevement of oil beads
        for i in range(number_of_beads):
            #Plot the trajectories for each bead i being bead and j being steps +2 to include all steps up to the last
            ax.plot(movement_x[:j+2, i], movement_y[:j+2, i], alpha=0.8, linewidth=0.5,zorder=1)
            
        
        #ax.set_title(f"Step number:{j+1}")
        ax.set_xlim(100,250)
        ax.set_ylim(100,250)
        ax.invert_yaxis()
        plt.show()


# Get minimum ice concentration for each bead across all timesteps
min_ice_per_bead = np.min(iceconc_array, axis=0)  # axis=0 = minimum over timesteps

for i in range(number_of_beads):
    print("Bead", i+1, "minimum ice concentration:" , round(min_ice_per_bead[i],2))