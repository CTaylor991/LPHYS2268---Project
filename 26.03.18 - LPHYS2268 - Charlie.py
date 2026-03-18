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

#2014-2023
start_yr = 2017

start_ensemble_number = 2024-start_yr

number_of_beads = 100



sea_ice_min = 0.15

# Central arctic
x0=175
y0=155

xlim_top = 260
xlim_bottom = 90

ylim_top = 260
ylim_bottom = 90



#Max number of files to cycke through
number_of_files = 3
file_range = range(number_of_files)
steps = 365 - start_date
total_steps = number_of_files * (365 - start_date)
#Required for tracking paths after each loop
long_running_step = 0

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
 
#initialise arrays to store movement rows are time steps and xolumns are beads
movement_x = np.zeros((total_steps + 1, number_of_beads))
movement_y = np.zeros((total_steps + 1, number_of_beads))

#initialise arrays to store ice concentration rows are time steps and xolumns are ice con
iceconc_array = np.zeros((total_steps + 1, number_of_beads))

#Store the initial position
movement_x[0] = bead_x
movement_y[0] = bead_y

#Array to capture those that drop out and when
drop_outs = np.zeros((number_of_beads,2))

#Array to track the active beads so we can skip the non-active ones in the movement
active_beads = np.zeros(number_of_beads)

#Plotting outside the loop rather than in it so it only happens once
 #----------Build custom colourmap-----------------------------
cmap = plt.cm.Blues_r
#land colours defined as nan
cmap.set_bad(color='darkgreen')   
#Set colour for when threshold is broken (open ocean, no sea ice)
cmap.set_under(color='midnightblue')

fig, ax = plt.subplots()
siconc_first = xr.open_dataset(f"siconc_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_yr}0601_{start_yr+1}0531_icemod.nc")['siconc']
siconc_first.isel(time_counter=start_date).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
ax.invert_yaxis()
ax.scatter(bead_x, bead_y, color='black', s=25)
ax.scatter(x0, y0, marker='*', color='sienna', s=100)
ax.set_title("Starting condition")
plt.show()


for year in file_range:
    #Start from November 1st in first file
    if year == 0:
        file_start_day = start_date 
    #start from beginning of file in any other file
    else:
        file_start_day = 0           # start from beginning of subsequent files
    
    steps = 365 - file_start_day
    
    #-------Loading in NetCDF files and plotting for background--------------------------------------
    start_year = year + start_yr
    finish_year = start_year + 1
    #pull the files, the :02d turns the value into two figures so we get 08, 09 and 10
    ds_i = xr.open_dataset(f"siconc_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
    ds_u = xr.open_dataset(f"sivelu_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
    ds_v = xr.open_dataset(f"sivelv_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
    siconc = ds_i['siconc']   
    
    #ensemble number decreases from 2014 to 2024
    start_ensemble_number = start_ensemble_number-1
    
    
    for i in range(number_of_beads):
        if active_beads[i] == 1:
            continue
        
        grid_x = int(bead_x[i])
        grid_y = int(bead_y[i])
        iceconc_array[0, i] = siconc.isel(time_counter=file_start_day, y=grid_y, x=grid_x).values
    
    
    for j in range(steps):
        #Break the loops if all beads have melted out or hit the land
        if np.all(active_beads == 1):
            print("All beads dropped out at step", long_running_step + j)
            break
        #Check if bead has melted out or hit land and then skip if it has
        for i in range(number_of_beads):
            if active_beads[i] == 1:
                continue
            
            init_x = bead_x[i]
            init_y = bead_y[i]
    
                    
            grid_x = int(init_x)
            grid_y = int(init_y)
            
            u_velocity = ds_u["sivelu"]
            u_velocity = u_velocity.isel(time_counter=file_start_day+j, y=grid_y, x=grid_x)
            
            v_velocity = ds_v["sivelv"]
            v_velocity = v_velocity.isel(time_counter=file_start_day+j, y=grid_y, x=grid_x)
            
            ice_conc = siconc.isel(time_counter=file_start_day+j, y=grid_y, x=grid_x)
            iceconc_array[long_running_step+j+1, i] = ice_conc.values
            
            #check if bead is either below ice conc value or hitting land, add to list and then skip movement if true
            if ice_conc < sea_ice_min or np.isnan(u_velocity.values) or np.isnan(v_velocity.values):
                active_beads[i] = 1
                if drop_outs[i, 0] == 0:
                    drop_outs[i, 0] = i+1
                    drop_outs[i, 1] = long_running_step + j + 1
                    #np.append(drop_outs, long_running_step, axis=1)
                    print("Bead", i+1, "dropped out at", long_running_step + j + 1)
                continue
                
            init_x = bead_x[i]
            init_y = bead_y[i]
            #Move along (86,400 seconds in a day, 25km grid size)
            bead_x[i] = bead_x[i] + (u_velocity.values * 86400/25000)
            bead_y[i] = bead_y[i] + (v_velocity.values * 86400/25000)
        
                   
            
            
        #Save the new position to the trajectory arrays. +1 because it is the next step in the array and to preserve x0 and y0
        movement_x[j+1+long_running_step] = bead_x
        movement_y[j+1+long_running_step] = bead_y
        
        
        if (j+1) % 20 == 0:
    
            fig, ax = plt.subplots()
            
            # Background sea ice for this timestep
            siconc.isel(time_counter=file_start_day+j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
            
            
            #Plot the new position of the beads and the centre point
            ax.scatter(x0,y0,marker='*',s=100,color='sienna',zorder=3)    
            ax.scatter(bead_x,bead_y,color='black',s=25,zorder=2)
            
        
            #Plotting of moevement of oil beads
            for i in range(number_of_beads):
                #Plot the trajectories for each bead i being bead and j being steps +2 to include all steps up to the last
                ax.plot(movement_x[:long_running_step+j+2, i], movement_y[:long_running_step+j+2, i], alpha=0.8, linewidth=0.5,zorder=1)
                
            
            #ax.set_title(f"Step number:{j+1}")
            ax.set_xlim(xlim_bottom,xlim_top)
            ax.set_ylim(ylim_bottom,ylim_top)
            
            ax.invert_yaxis()
            plt.show()
    
     #Move long_running step along
    long_running_step = long_running_step + steps

# 1 final plot
fig, ax = plt.subplots()

# Background sea ice for this timestep
siconc.isel(time_counter=file_start_day+j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)


#Plot the new position of the beads and the centre point
ax.scatter(x0,y0,marker='*',s=100,color='sienna',zorder=3)    
ax.scatter(bead_x,bead_y,color='black',s=25,zorder=2)


#Plotting of moevement of oil beads
for i in range(number_of_beads):
    #Plot the trajectories for each bead i being bead and j being steps +2 to include all steps up to the last
    ax.plot(movement_x[:long_running_step+j+2, i], movement_y[:long_running_step+j+2, i], alpha=0.8, linewidth=0.5,zorder=1)
    

#ax.set_title(f"Step number:{j+1}")
ax.set_xlim(xlim_bottom,xlim_top)
ax.set_ylim(ylim_bottom,ylim_top)

ax.invert_yaxis()
plt.show()

# Get minimum ice concentration for each bead across all timesteps
min_ice_per_bead = np.min(iceconc_array, axis=0)  # axis=0 = minimum over timesteps
print("List of drop outs\nBead, step\n", drop_outs)

# for i in range(number_of_beads):
#     print("Bead", i+1, "minimum ice concentration:" , round(min_ice_per_bead[i],2))