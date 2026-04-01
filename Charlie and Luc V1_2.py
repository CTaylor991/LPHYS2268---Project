# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 14:33:38 2026

@author: anegu
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import xarray as xr

#Import tool to calculate the convexhull around active beads
#To estimate the occupied surface area
from scipy.spatial import ConvexHull

#--------------Variables that may be changed----------------------------------
#Note 1st of November is j = 153
start_date = 153

#Year to start from (2014-2023)
start_yr = 2014


#Total number of beads
number_of_beads = 100
#Set threshold for sea ice conc
sea_ice_min = 0.15

# Central arctic coordinates
x0=175
y0=155

#Set limits for graph figures
xlim_top = 265
xlim_bottom = 85
ylim_top = 265
ylim_bottom = 85
#----------------------------------------------------------------------------

#--------------Varations for functionality----------------------------
#Align starting year to ensemble number 2014 = 10, 2015 = 09
start_ensemble_number = 2024-start_yr

#Max number of files to cycke through
number_of_files = 10
file_range = range(number_of_files)

#Number of steps for first year and for total length of programme
steps = 365 - start_date
total_steps = number_of_files * (365 - start_date)

#Required for tracking paths after each loop
long_running_step = 0
#Required for plotting final plot
last_active_step = 0
#-------------------------------------------------------------------------

#-----------------------Create bead formation------------------------
#Set side length
n = int(np.sqrt(number_of_beads))
#Create 10 evenly spaced data points that are centred on 0 and offset from the edges by 0.5
spacing = np.linspace(-0.5, 0.5, n, endpoint=False) + 0.5/n
    
#Create a 10x10 square grid shape
gx, gy = np.meshgrid(spacing, spacing)  
    
#shape the grid about the coordinate points
bead_x = x0 + gx.ravel()  
bead_y = y0 + gy.ravel() 
#------------------------------------------------------------------------- 

#----------------Arrays for storing data throughout the script-------------------
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


#------Build custom colourmap and plot the total area at day 0-----------------------------
cmap = plt.cm.Blues_r
#land colours defined as nan
cmap.set_bad(color='darkgreen')   
#Set colour for when threshold is broken (open ocean, no sea ice)
cmap.set_under(color='midnightblue')

#Plot inital map of total area with highlighted square of what we will look at
fig, ax = plt.subplots()
siconc_first = xr.open_dataset(f"siconc_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_yr}0601_{start_yr+1}0531_icemod.nc")['siconc']
siconc_first.isel(time_counter=start_date).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
#Invert y axis to make the plot look more familiar
ax.invert_yaxis()
#Plot highlighting square
ax.plot([xlim_bottom, xlim_top, xlim_top, xlim_bottom, xlim_bottom], [ylim_bottom, ylim_bottom, ylim_top, ylim_top, ylim_bottom], color='r')
#plot the beads (probably not needed)
ax.scatter(bead_x, bead_y, color='black', s=25)
#plot the origin of the beads
ax.scatter(x0, y0, marker='*', color='sienna', s=100)
ax.set_title("Starting condition")
plt.show()
#------------------------------------------------------------------------

#----------------IDW fuction-----------------------------------------------
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

#-------------------------------------------------------------------------------

#---------------Function to calculate occupied surface area------------------------------
def convex_hull_area_km2(x_positions, y_positions, cell_size_km=25):
    #Combine points coordinates into a list of 2D points
    points = np.column_stack((x_positions, y_positions))

    # 03 points are needed to form a surface
    if len(points) < 3:
        return 0.0

    try:
        #Convex hull around active beads
        hull = ConvexHull(points)
        area_grid_units = hull.volume
        return area_grid_units * (cell_size_km ** 2)
    except:
        return 0.0
    
#............Save surfaces....................................................
#Store daily occupied area for each ensemble and to list their label for each yr
all_areas = []
all_labels = []

#---------------Loop to run through to do data analysis and plotting--------------
for year in file_range:
    #Reset bead positions at the start of each ensemble
    bead_x = x0 + gx.ravel()
    bead_y = y0 + gy.ravel()
    #All beads are active to de new ensemble
    active_beads = np.zeros(number_of_beads)
   
    
    # List for local year
    area_this_year = []
    
    #Break the loops if all beads have melted out or hit the land
    if np.all(active_beads == 1):
        break
    
    #Start from November 1st in first file
    if year == 0:
        file_start_day = start_date 
    #start from beginning of file in any other file
    else:
        file_start_day = start_date           
    
    #Set step to = 0 for first day
    steps = 365 - file_start_day
    
    #-------Loading in NetCDF files and plotting for background--------------------------------------
    start_year = year + start_yr
    finish_year = start_year + 1
    #pull the files, the :02d turns the value into two figures so we get 08, 09 and 10
    ds_i = xr.open_dataset(f"siconc_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
    ds_u = xr.open_dataset(f"sivelu_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
    ds_v = xr.open_dataset(f"sivelv_sipn_easegrid_2024-06-01_{start_ensemble_number:02d}_1d_{start_year}0601_{finish_year}0531_icemod.nc")
    siconc = ds_i['siconc']   
    
    # Calculation of the initail surface of each year using all active beads
    active_x = bead_x[active_beads == 0]
    active_y = bead_y[active_beads == 0]
    
    #Store the occupied area at day 0
    area0 = convex_hull_area_km2(active_x, active_y, cell_size_km=25)
    area_this_year.append(area0)


    #ensemble number decreases from 2014 to 2024
    start_ensemble_number = start_ensemble_number-1
    

    #Store initial ice concentration using IDW
    for i in range(number_of_beads):
        #Check if bead is active and skip if not
        if active_beads[i] == 1:
            continue
        #store initial ice concentration value
        iceconc_array[long_running_step, i] = idw_interpolation(siconc, file_start_day, bead_x[i], bead_y[i])    
    
    #Run through the steps to do data processing
    for j in range(steps):
        #Break the loops if all beads have melted out or hit the land
        if np.all(active_beads == 1):
            print("All beads dropped out at step", long_running_step + j)
            break
        
        #Check if bead has melted out or hit land and then skip if it has
        for i in range(number_of_beads):
            if active_beads[i] == 1:
                continue
            
            #Calculate using IDW the two velocities and the sea ice concentration
            u_velocity = idw_interpolation(ds_u["sivelu"], file_start_day+j, bead_x[i], bead_y[i])
            v_velocity = idw_interpolation(ds_v["sivelv"], file_start_day+j, bead_x[i], bead_y[i])
            ice_conc = idw_interpolation(siconc, file_start_day+j, bead_x[i], bead_y[i])
            
            #add value to array
            iceconc_array[long_running_step+j+1, i] = ice_conc

            
            #check if bead is either below ice conc value or hitting land, add to list and then skip movement if true
            if ice_conc < sea_ice_min or np.isnan(u_velocity) or np.isnan(v_velocity):
                active_beads[i] = 1
                if drop_outs[i, 0] == 0:
                    drop_outs[i, 0] = i+1
                    drop_outs[i, 1] = long_running_step + j + 1
                    #np.append(drop_outs, long_running_step, axis=1)
                    print("Bead", i+1, "dropped out at", long_running_step + j + 1)
                continue
                
            
            #Move along (86,400 seconds in a day, 25km grid size)
            bead_x[i] = bead_x[i] + (u_velocity * 86400/25000)
            bead_y[i] = bead_y[i] + (v_velocity * 86400/25000)
        
                   
            
            
        #Save the new position to the trajectory arrays. +1 because it is the next step in the array and to preserve x0 and y0
        movement_x[j+1+long_running_step] = bead_x
        movement_y[j+1+long_running_step] = bead_y
        
        # Select only beads that are still active
        active_x = bead_x[active_beads == 0]
        active_y = bead_y[active_beads == 0]
        
        #Occupied are for the current day
        area_today = convex_hull_area_km2(active_x, active_y, cell_size_km=25)
        #Store the occupied area for this ensemble
        area_this_year.append(area_today)

        
        #keep moving along unless all beads are melted out or hit land
        if not np.all(active_beads == 1):
            last_active_step = long_running_step + j + 1
        
        #Plot to terminal depending on whatever interval we want
        if (j+1) % 60 == 0:
            print("Working on step:", j+1 ," to ", j+61)
            fig, ax = plt.subplots()
            
            #Background sea ice for this timestep
            siconc.isel(time_counter=file_start_day+j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)
            
            
            #Plot the new position of the beads and the centre point
            ax.scatter(x0,y0,marker='*',s=100,color='sienna',zorder=3)   
            
            #Split beads into active and inactive beads for plotting
            active_x = bead_x[active_beads == 0]
            active_y = bead_y[active_beads == 0]
            
            inactive_x = bead_x[active_beads == 1]
            inactive_y = bead_y[active_beads == 1]
            
            #Plot active beads as dots
            ax.scatter(active_x, active_y, color='black', s=25, zorder=2)
            
            #Plot inactive beads as crosses
            ax.scatter(inactive_x, inactive_y, color='red', marker='x', s=25, zorder=2)
                        
        
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
    
    # Keep the full area evolution curve for this ensemble and the labels of those ensembles
    all_areas.append(area_this_year)
    all_labels.append(f"{start_year}-{finish_year}")


#------------------------ 1 final plot--------------------
fig, ax = plt.subplots()

# Background sea ice for this timestep
siconc.isel(time_counter=file_start_day+j).plot(ax=ax, cmap=cmap, vmin=sea_ice_min, vmax=1)

#Plot the new position of the beads and the centre point
ax.scatter(x0,y0,marker='*',s=100,color='sienna',zorder=3)    
#Split beads into active and inactive beads for plotting
active_x = bead_x[active_beads == 0]
active_y = bead_y[active_beads == 0]

inactive_x = bead_x[active_beads == 1]
inactive_y = bead_y[active_beads == 1]

#Plot active beads as dots
ax.scatter(active_x, active_y, color='black', s=25, zorder=2)

#Plot inactive beads as crosses
ax.scatter(inactive_x, inactive_y, color='red', marker='x', s=25, zorder=2)

#Plotting of moevement of oil beads
for i in range(number_of_beads):
    #Plot the trajectories for each bead i being bead and j being steps +2 to include all steps up to the last
    ax.plot(movement_x[:last_active_step, i], movement_y[:last_active_step, i], alpha=0.8, linewidth=0.5, zorder=1)

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

#----------------Plot the occuppied area for all ensembles------------------------------
plt.figure(figsize=(10,7))

for k in range(len(all_areas)):
    days = np.arange(len(all_areas[k]))
    plt.plot(days, all_areas[k], label=all_labels[k])
    

#-----------Mean trend--------------------
min_length = min(len(a) for a in all_areas)
areas_array = np.array([a[:min_length] for a in all_areas])
mean_area = np.mean(areas_array, axis=0)

days = np.arange(min_length)
plt.plot(days, mean_area, color='black', linewidth=3, label='Mean')

plt.xlabel("Simulation day")
plt.ylabel("Occupied area (km²)")
plt.title("Variation of occupied area for each ensemble")
plt.grid(True)
plt.legend()
plt.show()