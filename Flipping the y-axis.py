import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("siconc_sipn_easegrid_2024-06-01_01_1d_20230601_20240531_icemod.nc")

sea_ice = ds["siconc"]

ax = sea_ice.isel(time_counter=20).plot()

# Flip the displayed y-axis only
ax.axes.invert_yaxis()

plt.show()