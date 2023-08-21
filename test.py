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

def extracting_SI_velocity(lat_range = [64,85], lon_range = [-40,20]):

    print("\n###########################################")
    print("##### - Extracting SI velocity data - #####")
    print("###########################################\n")

    u_si_file = "C:/Users/augustin/downloads/sea_ice_drift.txt"

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    u_ds = xr.open_dataset(u_si_file, decode_times = False)
    print(u_ds)
    usi= u_ds['usi'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude < lat_max))
    
    lon = v_ds['longitude']
    lat = v_ds['latitude']
    time =  v_ds['time']

    v_ds.close
    u_ds.close
    print("############ - Data extracted - ###########")
    return lon, lat, usi, vsi, time
extracting_SI_velocity()