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
    """  Return the date given a number n in time[n] in [0,132]
    time[0] = 596148
    time[-1] = 618660
    """
    return date(1950,1,1) + timedelta(hours=596148 + time_snapshot * 24)

def give_time_snapshot(year,month,day):
    """   
        Return the n indix of time[n] given a date. Bug to fix : give_date(give_time_snapshot(date)) doesn't return exactly the date.
    """
    days_from_start = (date(year,month,day) - date(2010,10,1)).days
    time_coverage = time[-1] - time[0] # In days
    time_step = time_coverage/len(time) # In days
    time_snapshot = days_from_start/time_step
    return int(time_snapshot)

def extracting_data(file = "C:/Users/Augustin/Downloads/Salinity_and_density", lat_range = [64,85], lon_range = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range

        time range: 03/01/2018 --> 07/29/2023
        output are DataArray.
    """
    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]
    ds = xr.open_dataset(file, decode_times = False)
    lon = ds['lon']
    lat = ds['lat']
    sal = ds['sos'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max))
    sal_uncertainty = ds['sos_error'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max))
    den = ds['dos'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max))
    den_uncertainty = ds['dos_error'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max))
    time =  ds['time']
    sal = sal.sel(depth=0)
    den = den.sel(depth=0)
    ds.close
    return lon, lat, sal, sal_uncertainty, den, den_uncertainty, time

def plot_mensual_salinity_mean(projection = ccrs.LambertConformal(central_longitude = -20), figsize = (9,7)):
    """
        Save plot of density and salinity map for march and sept.
    """
    dates = [[2018,3,15],[2018,9,15],[2019,3,15]]
    for date in dates:
        print(f"### - Saving salinity field: {date[0]}-{date[1]} - ###\n")
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
        axs.set_extent([-70, 25, 55, 85], crs = ccrs.PlateCarree())
        axs.coastlines()
        axs.gridlines()
        levels = np.linspace(29.5,35.2,1000)
        cs = axs.contourf(lon, lat, mensual_salinity_mean(date[0],date[1]), levels = levels, cmap = "cmo.haline", transform=ccrs.PlateCarree())
        axs.set_title("Mensual mean Sanility value in psu for {}/{}".format(date[0],date[1]))
        fig.colorbar(cs, ax = axs)
        plt.savefig(f"Plots/mean/Salinity/salinity_mean_{date[0]}-{date[1]}.png")

def mensual_salinity_mean(year, month):
    """
        Return the mensual mean value of den and sal
    """
    useful_date = []
    for time_ in time:
        corresponding_date = date(1950,1,1) + timedelta(hours = int(time_))
        if corresponding_date.year == year and corresponding_date.month == month:
            useful_date.append(int(time_))
            print(corresponding_date)
    recorded_sal = np.array([sal.sel(dict(time=n)) for n in useful_date])
    mean_sal = np.nanmean(recorded_sal, axis = 0)
    return mean_sal

def mensual_density_mean(year, month):
    """
        Return the mensual mean value of den and sal
    """
    useful_date = []
    for time_ in time:
        corresponding_date = date(1950,1,1) + timedelta(hours = int(time_))
        if corresponding_date.year == year and corresponding_date.month == month:
            useful_date.append(int(time_))
            print(corresponding_date)
    recorded_den = np.array([den.sel(dict(time=n)) for n in useful_date])
    mean_den = np.nanmean(recorded_den, axis = 0)
    return mean_den

def plot_mensual_density_mean(projection = ccrs.LambertConformal(central_longitude = -20), figsize = (9,7)):
    """
        Save plot of density and salinity map for march and sept.
    """
    dates = [[2018,3,15],[2018,9,15],[2019,3,15]]
    for date in dates:
        print(f"### - Saving density field: {date[0]}-{date[1]} - ###\n")
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
        axs.set_extent([-70, 25, 55, 85], crs = ccrs.PlateCarree())
        axs.coastlines()
        axs.gridlines()
        levels = np.linspace(1023.5,1028.25,1000)
        cs = axs.contourf(lon, lat, mensual_density_mean(date[0],date[1]), levels = levels, cmap = "cmo.dense", transform=ccrs.PlateCarree())
        axs.set_title("Mensual mean density value in kg/m^3 for {}/{}".format(date[0],date[1]))
        fig.colorbar(cs, ax = axs)
        plt.savefig(f"Plots/mean/Density/density_mean_{date[0]}-{date[1]}.png")
lon, lat, sal, sal_uncertainty, den, den_uncertainty, time = extracting_data()
plot_mensual_density_mean()
