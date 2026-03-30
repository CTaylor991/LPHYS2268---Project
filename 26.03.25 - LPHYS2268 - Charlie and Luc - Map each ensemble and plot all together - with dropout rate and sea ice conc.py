# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:14:20 2026

Notes: grid is 25km x 25km, it is 361 x 361 grid cells.

@author: charl
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import xarray as xr
import imageio
import os

#Variables------------------------------------------
#Note 1st of November is j = 153
start_date = 153

#2014-2023
#start_yr = 2017
Year_range = [2014,2015,2016,2017,2018,2019,2020,2021,2022,2023]
#Year_range = [2014,2015]

#Dictionaries for final plots
x_movement = {}
y_movement = {} 
drop_outs = {}
last_active_steps = {}
sea_ice_concs = {}

#Data path
DATA_PATH = r"/home/elic/ctaylor/LPHYS2268/Data_files"

# IDW fuction
def idw_interpolation(var, time_index, x_pos, y_pos):

    x0 = int(np.floor(x_pos))
    x1 = x0 + 1
    y0 = int(np.floor(y_pos))
    y1 = y0 + 1

    points = [(x0,y0),(x1,y0),(x0,y1),(x1,y1)]

    weighted_sum = 0
    weight_sum = 0

    for px,py in points:

        value = var.isel(time_counter=time_index,y=py,x=px).values

        dist = np.sqrt((x_pos-px)**2 + (y_pos-py)**2)

        if dist == 0:
            return value

        weight = 1/(dist**2)

        weighted_sum += weight*value
        weight_sum += weight

    return weighted_sum/weight_sum

for start_yr in Year_range:
    constant_start_yr = start_yr
    start_ensemble_number = 2024-start_yr
    print(f"Working on start year: {constant_start_yr}")
    
    number_of_beads = 100
      
    sea_ice_min = 0.15
    
    # Central arctic
    x0=175
    y0=155
    
    xlim_top = 270
    xlim_bottom = 90
    
    ylim_top = 280
    ylim_bottom = 90
    
    
    #Max number of files to cycke through
    number_of_files = 4
    file_range = range(number_of_files)
    steps = 365 - start_date
    total_steps = number_of_files * (365 - start_date)
    #Required for tracking paths after each loop
    long_running_step = 0
    #Required for plotting final plot
    last_active_step = 0

    ice_area_list = []


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
    drop_outs_array = np.zeros((number_of_beads,2))
    
    #Array to track the active beads so we can skip the non-active ones in the movement
    active_beads = np.zeros(number_of_beads)
    
    #Plotting outside the loop rather than in it so it only happens once
     #----------Build custom colourmap-----------------------------
    cmap = plt.cm.Blues_r
    #land colours defined as nan
    cmap.set_bad(color='darkgreen')   
    #Set colour for when threshold is broken (open ocean, no sea ice)
    cmap.set_under(color='midnightblue')
        
    
        
    for year in file_range:
        #Break the loops if all beads have melted out or hit the land
        if np.all(active_beads == 1):
            break
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
        #have to check the file because in the 2023 file it tries to run into 2024 which we dontt have data for and it breaks the code, this way it will just stop the loop when it gets to that point
        file_i = f"{DATA_PATH}/siconc_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc"
        file_u = f"{DATA_PATH}/sivelu_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc"
        file_v = f"{DATA_PATH}/sivelv_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc"
        #This does the actual chck and then breaks the loop if needed
        if not os.path.exists(file_i) or not os.path.exists(file_u) or not os.path.exists(file_v):
            print(f"File not found for year {start_year}-{finish_year}, stopping loop.")
            break
        
        ds_i = xr.open_dataset(f"{DATA_PATH}/siconc_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
        ds_u = xr.open_dataset(f"{DATA_PATH}/sivelu_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
        ds_v = xr.open_dataset(f"{DATA_PATH}/sivelv_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
        siconc = ds_i['siconc']   
        print(siconc.shape)
        
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
        
                u_velocity = idw_interpolation(ds_u["sivelu"], file_start_day+j, bead_x[i], bead_y[i])
                v_velocity = idw_interpolation(ds_v["sivelv"], file_start_day+j, bead_x[i], bead_y[i])
                ice_conc = idw_interpolation(siconc, file_start_day+j, bead_x[i], bead_y[i])
                
                iceconc_array[long_running_step+j+1, i] = ice_conc
    
                #check if bead is either below ice conc value or hitting land, add to list and then skip movement if true
                if ice_conc < sea_ice_min or np.isnan(u_velocity) or np.isnan(v_velocity):
                    active_beads[i] = 1
                    if drop_outs_array[i, 0] == 0:
                        drop_outs_array[i, 0] = i+1
                        drop_outs_array[i, 1] = long_running_step + j
                        #print("Bead", i+1, "dropped out at", long_running_step + j + 1)
                    continue
                    
                init_x = bead_x[i]
                init_y = bead_y[i]
                #Move along (86,400 seconds in a day, 25km grid size)
                bead_x[i] = bead_x[i] + (u_velocity * 86400/25000)
                bead_y[i] = bead_y[i] + (v_velocity * 86400/25000)
            #gris is 25kmx25km so multiply the number of grid cells with ice by 625 to get area in km^2
            ice_area_list.append(float(siconc.isel(time_counter=file_start_day+j).sum()) * 625)
            #Save the new position to the trajectory arrays. +1 because it is the next step in the array and to preserve x0 and y0
            movement_x[j+1+long_running_step] = bead_x
            movement_y[j+1+long_running_step] = bead_y
            

            
            # x_movement_constant_start_yr = []
            # y_movement_constant_start_yr = [] 
            
            if not np.all(active_beads == 1):
                last_active_step = long_running_step + j + 1
            
          
         #Move long_running step along
        long_running_step = long_running_step + steps
    
        #add values to dictionary
    #print(f"Year {constant_start_yr}: last_active_step = {last_active_step}, long_running_step = {long_running_step}")
    x_movement[constant_start_yr] = movement_x[:last_active_step, :].copy()
    y_movement[constant_start_yr] = movement_y[:last_active_step, :].copy()
    drop_outs[constant_start_yr] = drop_outs_array.copy()
    last_active_steps[constant_start_yr] = last_active_step
    #Add the sea ice concentration for this step to the array
    sea_ice_concs[constant_start_yr] = ice_area_list.copy()
    #------------------------ 1 final plot--------------------
    fig, ax = plt.subplots()
    siconc.isel(time_counter=file_start_day+j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
    ax.scatter(x0,y0,marker='*',s=100,color='k',edgecolors='white',zorder=3)    
    active_x = bead_x[active_beads == 0]
    active_y = bead_y[active_beads == 0]
    inactive_x = bead_x[active_beads == 1]
    inactive_y = bead_y[active_beads == 1]
    ax.scatter(active_x, active_y, color='black', s=25, zorder=2)
    ax.scatter(inactive_x, inactive_y, color='red', marker='x', s=25, zorder=4)
    for i in range(number_of_beads):
        ax.plot(movement_x[:last_active_step, i], movement_y[:last_active_step, i], alpha=0.8, linewidth=0.5, zorder=1)
    ax.set_xlim(xlim_bottom,xlim_top)
    ax.set_ylim(ylim_bottom,ylim_top)
    ax.set_title(f"Final state for start year: {constant_start_yr}, after {last_active_step} days, start point ({x0}, {y0})")
    ax.invert_yaxis()
    fig.savefig(f"Final_state_{constant_start_yr} , start point ({x0}, {y0}).jpg", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Get minimum ice concentration for each bead across all timesteps
    min_ice_per_bead = np.min(iceconc_array, axis=0)
    print("List of drop outs\nBead, step\n", drop_outs)

#Map plot of entire ensemble+
fig, ax = plt.subplots()
#load in siconc file (doesn't matter what year because we are going to remove ice values and land shouldn't move)
siconc_background = xr.open_dataset(f"{DATA_PATH}/siconc_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")['siconc'].isel(time_counter=0)

#set all ice values to 0
masked = siconc_background * 0
masked.plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)

#colour range for each ensemble run
colours = plt.cm.tab10(np.linspace(0, 1, len(Year_range)))

#run through and plot the ensemble members
for idx, yr in enumerate(Year_range):
    for i in range(number_of_beads):
        ax.plot(x_movement[yr][:, i], y_movement[yr][:, i],alpha=0.5, linewidth=0.5, color=colours[idx])
        dropout_step = int(drop_outs[yr][i, 1])
        if dropout_step > 0: 
            idx_step = min(dropout_step, last_active_steps[yr] - 1)
            ax.scatter(x_movement[yr][idx_step, i], y_movement[yr][idx_step, i],color='k', marker='x', edgecolors='white', s=25, zorder=4)
    ax.plot([], [], color=colours[idx], label=f"{yr} ({last_active_steps[yr]} days)")
        
ax.scatter(x0, y0, marker='*', color='k',edgecolors='white', s=100, zorder=3)
ax.set_xlim(xlim_bottom, xlim_top)
ax.set_ylim(ylim_bottom, ylim_top)
ax.invert_yaxis()
ax.legend(title="Start year")
ax.set_title(f"Ensemble of trajectories across all starting years, start point ({x0}, {y0})")
fig.savefig(f"Ensemble_plot_all_years, start point ({x0}, {y0}).jpg", dpi=300, bbox_inches='tight')
plt.show()

#Dropout plot for entire ensemble
fig, ax = plt.subplots()
for idx, yr in enumerate(Year_range):
    
    drop_outs_year = drop_outs[yr]
    dropout_days = drop_outs_year[:, 1]
    # remove beads that never dropped out or they will be plotted at day 0 and skew the plot
    dropout_days = dropout_days[dropout_days > 0]  

    max_day = last_active_steps[yr]
    counts = np.zeros(max_day + 2)  # +2 so that the dropout graph goes back to zero after the last dropout day

    for day in dropout_days:
        counts[int(day)] += 1

    ax.plot(range(max_day + 2), counts, color=colours[idx], label=str(yr))
ax.set_xlabel("Day")
ax.set_ylabel("Number of beads dropping out")
ax.set_title(f"Bead dropouts per day across ensemble members, start point ({x0}, {y0})")
ax.legend(title="Start year")
fig.savefig(f"Ensemble_dropouts_all_years, start point ({x0}, {y0}).jpg", dpi=300, bbox_inches='tight')
plt.show()

#-------------------------Plot of sea ice concentration across ensemble members----------------------
fig, ax = plt.subplots()
for idx, yr in enumerate(Year_range):
    ax.plot(sea_ice_concs[yr], color=colours[idx], label=str(yr))
ax.set_xlabel("Day")
ax.set_ylabel("Ice area (km²)")
ax.set_title(f"Sea ice area over time, start point ({x0}, {y0})")
ax.legend(title="Start year")
fig.savefig(f"Ensemble_sea_ice_concentration_all_years, start point ({x0}, {y0}).jpg", dpi=300, bbox_inches='tight')
plt.show()