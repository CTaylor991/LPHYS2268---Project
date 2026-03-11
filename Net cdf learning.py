# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:16:08 2026

@author: charl
"""

import xarray as xr
import matplotlib.pyplot as plt

# Open file
ds = xr.open_dataset("siconc_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")

# Print dataset summary
print("\nDataset summary--------------------------------------------------------")
print(ds)
print("--------------------------------------------------------")
# Dimensions
print("\n\n\nDimensions:--------------------------------------------------------")
print(ds.dims)
print("--------------------------------------------------------")
# Coordinates
print("\n\n\nCoordinates:--------------------------------------------------------")
print(ds.coords)
print("--------------------------------------------------------")
# Variables
print("\n\n\nVariables:--------------------------------------------------------")
print(ds.data_vars)
print("--------------------------------------------------------")
# Variable keys
print("\n\n\nVariable keys:--------------------------------------------------------")
print(ds.variables.keys())
print("--------------------------------------------------------")


print("\n\n\nSea Ice info:--------------------------------------------------------")
sea_ice = ds["siconc"]
print(sea_ice)
print("--------------------------------------------------------")


print("\n\n\nMetadata for dataset:--------------------------------------------------------")
print(ds.attrs)
print("--------------------------------------------------------")

print("\n\n\nMetadata for dataset:--------------------------------------------------------")
print(sea_ice.attrs)
print("--------------------------------------------------------")


#Plot sea ice at a set time as a heat map-------------------------------------

ax = sea_ice.isel(time_counter=220).plot()
ax.axes.invert_yaxis()

plt.show()

print(ds.time_counter.min().values)
print(ds.time_counter.max().values)

print(sea_ice.min().values)
print(sea_ice.max().values)


#Plot the sea ice fraction over the whole time period
point = sea_ice.isel(y=180, x=180)

point.plot()

plt.show()

#Mean across the whole time period for each individual spot plotted as a heat map
mean_ice = sea_ice.mean(dim="time_counter")

mean_ice.plot()

plt.show()

#Sea ice 
y_pos = 180
x_pos = 180
time = 200
value = sea_ice.isel(time_counter=time, y=y_pos, x=x_pos)

print("Ice percentage at position:",value.values)