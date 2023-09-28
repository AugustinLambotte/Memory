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

def extracting_data(area_lat = [75,77.5], area_lon = [-10,0]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range
    """
    print("\n##############################")
    print("##### - Extracting data - ####")
    print("##############################\n")

    lon_min = area_lon[0]
    lon_max = area_lon[1]
    lat_min = area_lat[0]
    lat_max = area_lat[1]

    v_file = "C:/Users/Augustin/Downloads/Northward_sea_water_velocity"
    u_file = "C:/Users/Augustin/Downloads/Eastward_sea_water_velocity"
    
    v_ds = xr.open_dataset(v_file, decode_times = False)
    u_ds = xr.open_dataset(u_file, decode_times = False)
    
    u_gos= u_ds['uo'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude < lat_max),drop = True)
    v_gos= v_ds['vo'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude < lat_max),drop = True)
    
    u_gos = u_gos.sel(depth = 0.494025)
    v_gos = v_gos.sel(depth = 0.494025)
    
    lon = v_ds['longitude']
    lat = v_ds['latitude']
    time =  v_ds['time']
    v_ds.close
    u_ds.close

    print('##### - Data extracted - #####\n')
    return lon, lat, u_gos, v_gos, time

def mensual_mean(year, month,u_gos,v_gos,time):

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

def index_cell(name ='Cell_A',area_lat = [75,77.5], area_lon = [-10,0]):

    lon, lat, u_gos,v_gos,time = extracting_data(area_lat=area_lat,area_lon=area_lon)
    kinetic_energy = []

    for year in range(2011,2020):
        kinetic_energy.append([])
        for month in range(1,13):
            print(f"### - computing kin ene cell {name[-1]}: {year}-{month} - ###\n")
            u, v = mensual_mean(year, month,u_gos,v_gos,time)
            
            kinetic_energy[-1].append(np.nanmean(1/2 * (u**2 + v**2))) #In [J/kg]
    kinetic_energy = np.array(kinetic_energy)
    np.savetxt("Data/"+name+"/Cell_"+name[-1]+"_index_monthly.txt",kinetic_energy)      
    
index_cell(name = 'Cell_A',area_lat=[77.5,80],area_lon = [-10,0])
index_cell(name = 'Cell_B',area_lat=[75,77.5],area_lon = [-10,0])
index_cell(name = 'Cell_C',area_lat=[72.5,75],area_lon = [-20,-10])
index_cell(name = 'Cell_D',area_lat=[70,72.5],area_lon = [-20,-10])
