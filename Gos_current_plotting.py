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
def give_date(time_snapshot): 
    """  Return the date given a number n in time[n] in [0,3936]
    time[0] = 532524
    time[-1] = 618708
    """
    return date(1950,1,1) + timedelta(hours=532524 + time_snapshot * 24)

def give_time_snapshot(year,month,day):
    """   Return the n indix of time[n] given a date. Bug to fix : give_date(give_time_snapshot(date)) doesn't return exactly the date.
    """
    days_from_start = (date(year,month,day) - date(2010,10,1)).days
    time_coverage = time[-1] - time[0] # In days
    time_step = time_coverage/len(time) # In days
    time_snapshot = days_from_start/time_step
    return int(time_snapshot)

def extracting_data(lat_range = [64,85], lon_range = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range
    """
    print("\n##############################")
    print("##### - Extracting data - ####")
    print("##############################\n")

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    v_file = "C:/Users/Augustin/Downloads/Northward_sea_water_velocity"
    u_file = "C:/Users/Augustin/Downloads/Eastward_sea_water_velocity"
    
    v_ds = xr.open_dataset(v_file, decode_times = False)
    u_ds = xr.open_dataset(u_file, decode_times = False)
    
    u_gos= u_ds['uo'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude < lat_max))
    v_gos= v_ds['vo'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude < lat_max))

    u_gos = u_gos.sel(depth = 0.494025)
    v_gos = v_gos.sel(depth = 0.494025)

    lon = v_ds['longitude']
    lat = v_ds['latitude']
    time =  v_ds['time']
    v_ds.close
    u_ds.close

    print('##### - Data extracted - #####\n')
    return lon, lat, u_gos, v_gos, time

def plot_gos_current_map_single_date(year, month, day, projection, figsize = (14,10)):
    """
        Plot a single plot of the sea_ice_thickness at the given date value.
    """
    time_snapshot = give_time_snapshot(year,month,day)
    fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize=figsize, subplot_kw={'projection': projection})
    current = [u_gos, v_gos]
    for ax, u_or_v in zip(axs.flat, current):
        ax.set_extent([-40, 20, 60, 85], crs = ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines()
        cs = ax.contourf(lon, lat, u_or_v[time_snapshot], transform=ccrs.PlateCarree())
        fig.colorbar(cs, ax = ax)

    axs[0].set_title("{}-{}. zonal gos current [m/s].".format(year, month))
    axs[1].set_title("{}-{}. meridional gos current [m/s].".format(year, month))        
    
    plt.savefig("Plots/Gos_current_{}-{}.png".format(year,month))

def plot_current(year,month,day,projection = ccrs.LambertConformal(central_longitude = -20), figsize = (14,10)):
    time_snapshot = give_time_snapshot(year,month,day)
    current_magnitude = np.sqrt(u_gos**2 + v_gos**2)
    fig,axs = plt.subplots(figsize = figsize, subplot_kw={'projection': projection})
    axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
    axs.coastlines()

    #Magnitude plot
    levels = np.linspace(0,0.56,10)
    cs = axs.contourf(lon,lat,current_magnitude[time_snapshot],cmap = "cmo.speed",levels = levels,transform =ccrs.PlateCarree())
    
    #Vector plot
    #mymap = plt.streamplot(np.array(lon[:]),np.array(lat[:]),np.array(u_gos[time_snapshot]),np.array(v_gos[time_snapshot]), transform=ccrs.PlateCarree(), density=4)
    axs.quiver(np.array(lon[:]),np.array(lat[:]),np.array(u_gos[time_snapshot]),np.array(v_gos[time_snapshot]), transform = ccrs.PlateCarree())
    fig.colorbar(cs, ax=axs)
    plt.show()

def save_plot_mar_sept(projection, figsize = (14,10), mean = True):
    dates = [[2011,3,15],[2011,9,15],[2012,3,15],[2012,9,15],[2013,3,15],[2013,9,15],[2014,3,15],[2014,9,15],[2015,3,15],[2015,9,15],[2016,3,15]
            ,[2016,9,15],[2017,3,15],[2017,9,15],[2018,3,15],[2018,9,15],[2019,3,15],[2019,9,15]]
    print("##################################\n")
    
    for date in dates:
        print(f"### - Saving gos: {date[0]}-{date[1]} - ###\n")
        time_snapshot = give_time_snapshot(date[0],date[1],date[2])        
        if mean == True:
            u_gos, v_gos = mensual_mean(date[0], date[1])
        else:
            u_gos, v_gos = u_gos[time_snapshot], v_gos[time_snapshot]
        
        current_magnitude = np.sqrt(u_gos**2 + v_gos**2)
        fig,axs = plt.subplots(figsize = figsize, subplot_kw={'projection': projection})
        axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
        axs.coastlines()
        axs.set_title(f"{date[0]}-{date[1]} Magnitude and streamplot of surface gos current [m/s]")
        #Magnitude plot
        levels = np.linspace(0,0.56,10)
        cs = axs.contourf(lon,lat,current_magnitude, cmap = "cmo.speed", levels = levels,transform =ccrs.PlateCarree())
    
        #Vector plot
        #mymap = plt.streamplot(np.array(lon[:]),np.array(lat[:]),np.array(u_gos),np.array(v_gos), transform=ccrs.PlateCarree(), density=4)
        
        skip = (slice(None,None,4),slice(None,None,4))
        axs.quiver(np.array(lon.isel(dict(longitude=slice(None,None,4)))),np.array(lat.isel(dict(latitude=slice(None,None,4)))),u_gos[skip],v_gos[skip], transform = ccrs.PlateCarree())
        fig.colorbar(cs, ax=axs)
        if mean == True:
            plt.savefig(f"Plots/mean/Gos_current/stream_plot_{date[0]}-{date[1]}.png")
        else:
            plt.savefig(f"Plots/Gos_current/stream_plot_{date[0]}-{date[1]}.png")
    print("############################################")



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
            print(corresponding_date)
        i+=1
    # Creation of an array with all the data during the month of interest
    recorded_ugos = np.array([u_gos.isel(time = n) for n in useful_index])
    recorded_vgos = np.array([v_gos.isel(time = n) for n in useful_index])
    mean_ugos = np.nanmean(recorded_ugos, axis = 0)
    mean_vgos = np.nanmean(recorded_vgos, axis = 0)

    return mean_ugos, mean_vgos



#### - Main - ####
lon, lat, u_gos,v_gos,time = extracting_data()
projection = ccrs.LambertConformal(central_longitude = -20)
save_plot_mar_sept(projection)
