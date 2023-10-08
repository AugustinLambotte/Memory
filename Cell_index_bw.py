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


"""
    Same as Cell_index_daily.py but we extract only data when there is cryosat sit data
"""
def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [60,83], lon_range = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range

        output are in DataArray.
    """

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]
    ds = xr.open_dataset(file, decode_times = False)
    lon = ds['Longitude'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    lat = ds['Latitude'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sit = ds['Sea_Ice_Thickness']#.where((ds.Sea_Ice_Thickness != 0)) 
    sit = sit.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sic = ds['Sea_Ice_Concentration']#.where((ds.Sea_Ice_Thickness != 0)) 
    sic = sic.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    time =  ds['Time']

    
    gate_marker = 0
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time, gate_marker

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
    print(u_ds.attrs)
    v_ds.close
    u_ds.close

    print('##### - Data extracted - #####\n')
    return lon, lat, u_gos, v_gos, time


def index_cell(name ='Cell_A',area_lat = [75,77.5], area_lon = [-10,0]):
    
    lon, lat, u_gos,v_gos,time = extracting_data(area_lat=area_lat,area_lon=area_lon)
    lon_sit,lat_sit, sit, sic, sit_uncertainty, sit_time, gate_marker= extracting_data_sit()

    starting_date = 734419
    rec_date_sit = []
    for time_ in sit_time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
        if corresponding_date.year >= year_ and corresponding_date.year < year_end:
            rec_date_sit.append(corresponding_date)


    # Retrieve all the time index corresponding to the month of interest.
    useful_index = []
    i = 0
    
    for time_ in time:
        corresponding_date = date(1950,1,1) + timedelta(hours = int(time_)) 
        if corresponding_date in rec_date_sit:
            useful_index.append(i)
        i+=1

    # Creation of an array with all the data during the year of interest
    recorded_ugos = np.array([u_gos.isel(time = n) for n in useful_index])
    recorded_vgos = np.array([v_gos.isel(time = n) for n in useful_index])
    kinetic_energy = [np.nanmean(1/2 * (recorded_ugos[day]**2 + recorded_vgos[day]**2)) for day in range(len(recorded_ugos))] #In [J/kg]
    print(len(kinetic_energy))

    np.savetxt("Data/bw/"+name+"/Cell_"+name[-1]+"_index_daily.txt",kinetic_energy)      

year_,year_end = 2011,2021
index_cell(name = 'Cell_A',area_lat=[77.5,80],area_lon = [-10,0])
index_cell(name = 'Cell_B',area_lat=[75,77.5],area_lon = [-10,0])
index_cell(name = 'Cell_C',area_lat=[72.5,75],area_lon = [-20,-10])
index_cell(name = 'Cell_D',area_lat=[70,72.5],area_lon = [-20,-10])
