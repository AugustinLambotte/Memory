from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from datetime import datetime, date, timedelta


####### - Functions definitions - #####


def give_date(time_snapshot): 
    """  
        Return the date given a number n in time[n] in [0,235]
        time_coverage_start = "2010-10-01"
        time_coverage_end = "2020-07-31" 
    """
    time_coverage = time[-1] - time[0] #In days
    time_step = time_coverage/len(time)
    days_from_start = time_step * time_snapshot

    return date(2010,10,1) + timedelta(days=int(days_from_start))

def give_time_snapshot(year,month,day,time):
    """   
        Return the n indix of time[n] given a date. Bug to fix : give_date(give_time_snapshot(date)) doesn't return exactly the date.
    """
    days_from_start = (date(year,month,day) - date(2010,10,1)).days
    time_coverage = time[-1] - time[0] # In days
    time_step = time_coverage/len(time) # In days
    time_snapshot = days_from_start/time_step
    return int(time_snapshot)

def extracting_SIT(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [64,85], lon_range = [-40,20]):
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
    lon = ds['Longitude']
    lat = ds['Latitude']
    sit = ds['Sea_Ice_Thickness'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))
    sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))
    time =  ds['Time']
    ds.close
    return lon, lat, sit, sit_uncertainty, time

def extracting_SI_velocity(lat_range = [64,85], lon_range = [-40,20]):

    print("\n###########################################")
    print("##### - Extracting SI velocity data - #####")
    print("###########################################\n")

    u_si_file = "C:/Users/augustin/downloads/eastward_sea_ice_velocity"
    v_si_file = "C:/Users/augustin/downloads/Northward_sea_ice_velocity"

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    u_ds = xr.open_dataset(u_si_file, decode_times = False)
    v_ds = xr.open_dataset(v_si_file, decode_times = False)
    usi= u_ds['usi'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude < lat_max))
    vsi= v_ds['vsi'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude < lat_max))

    lon = v_ds['longitude']
    lat = v_ds['latitude']
    time =  v_ds['time']

    v_ds.close
    u_ds.close
    print("############ - Data extracted - ###########")
    return lon, lat, usi, vsi, time

def extracting_GOS(lat_range = [64,85], lon_range = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, u_gos, v_gos and time 
        restricted on the area defined by lat_range and lon_range
    """
    print("\n##################################")
    print("##### - Extracting GOS data - ####")
    print("##################################\n")

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    v_file = "C:/Users/Augustin/Downloads/northward_cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_1690898821778.nc"
    u_file = "C:/Users/Augustin/Downloads/eastward_cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_1690898632385.nc"
    
    v_ds = xr.open_dataset(v_file, decode_times = False)
    u_ds = xr.open_dataset(u_file, decode_times = False)
    
    u_gos= u_ds['ugos'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude < lat_max))
    v_gos= v_ds['vgos'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude < lat_max))

    lon = v_ds['longitude']
    lat = v_ds['latitude']
    time =  v_ds['time']

    v_ds.close
    u_ds.close
    print('##### - Data extracted - #####\n')
    return lon, lat, u_gos, v_gos, time

def current_map(lon, lat, u_gos, v_gos, time, year, month, day,projection = ccrs.LambertConformal(central_longitude = -20), figsize = (14,10)):
    """
        Plot GOS current map given u_gos, v_gos, time, year, month and day
    """
    time_snapshot = give_time_snapshot(year,month,day,time)
    current_magnitude = np.sqrt(u_gos**2 + v_gos**2)
    fig,axs = plt.subplots(figsize = figsize, subplot_kw={'projection': projection})
    axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
    axs.coastlines()

    #Magnitude plot
    cs = axs.contourf(lon,lat,current_magnitude[time_snapshot],transform =ccrs.PlateCarree())
    
    #Vector plot
    mymap = plt.streamplot(np.array(lon[:]),np.array(lat[:]),np.array(u_gos[time_snapshot]),np.array(v_gos[time_snapshot]), transform=ccrs.PlateCarree(), density=4)
    #axs.quiver(np.array(lon[:]),np.array(lat[:]),np.array(u_gos[time_snapshot]),np.array(v_gos[time_snapshot]), transform = ccrs.PlateCarree())
    fig.colorbar(cs, ax=axs)
    plt.show()

