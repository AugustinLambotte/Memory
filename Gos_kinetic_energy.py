from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from datetime import datetime, date, timedelta
import cmocean
import matplotlib.path as mpath
from Gos_current_plotting import extracting_data

def mensual_mean(year, month):

    # Retrieve all the time index corresponding to the month of interest.
    useful_date = []
    useful_index = []
    i = 0
    for time_ in time:
        corresponding_date = date(1950,1,1) + timedelta(hours = int(time_)) 
        if corresponding_date.year == year and corresponding_date.month == month:
            useful_date.append(int(time_))
            useful_index.append(i)
        i+=1
    # Creation of an array with all the data during the month of interest
    recorded_ugos = np.array([u_gos.isel(time = n) for n in useful_index])
    recorded_vgos = np.array([v_gos.isel(time = n) for n in useful_index])
    mean_ugos = np.nanmean(recorded_ugos, axis = 0)
    mean_vgos = np.nanmean(recorded_vgos, axis = 0)
    return mean_ugos, mean_vgos

lon, lat, u_gos,v_gos,time = extracting_data()
for year in range(2011,2021):
    for month in range(1,13):
        print(f"### - Saving gos: {year}-{month} - ###\n")
        u, v = mensual_mean(year, month)
        
        kinetic_energy = 1/2 * (u**2 + v**2) #In [J/kg]
        
        fig,axs = plt.subplots(figsize = (14,10), subplot_kw={'projection': ccrs.LambertConformal(central_longitude = -18)})
        #axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
        xlim = [-43, 16]
        ylim = [61, 81]
        lower_space = 3 
        rect = mpath.Path([[xlim[0], ylim[0]],
                        [xlim[1], ylim[0]],
                        [xlim[1], ylim[1]],
                        [xlim[0], ylim[1]],
                        [xlim[0], ylim[0]],
                        ]).interpolated(20)
        proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
        rect_in_target = proj_to_data.transform_path(rect)
        axs.set_boundary(rect_in_target)
        axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
        axs.coastlines()
        axs.gridlines()
        axs.set_title(f"{year}-{month} Surface gos current - Kinetic energy [J/kg]")
        #Magnitude plot
        levels = np.linspace(0,0.2,10)
        cs = axs.contourf(lon,lat,kinetic_energy, cmap = "cmo.speed",levels = levels,transform =ccrs.PlateCarree())
        #Vector plot
        #mymap = plt.streamplot(np.array(lon[:]),np.array(lat[:]),np.array(u_gos),np.array(v_gos), transform=ccrs.PlateCarree(), density=4)
        fig.colorbar(cs, ax=axs, ticks = [0,0.1,0.2,0.3,0.4])
        plt.savefig(f"Plots/mean/Kinetic_energy/{year}/GOS_kin_{year}-{month}.png")
