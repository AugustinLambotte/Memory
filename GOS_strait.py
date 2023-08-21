from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from datetime import datetime, date, timedelta

from Useful_function import extracting_GOS, give_date, give_time_snapshot, current_map

lon,lat,u_gos,v_gos,time = extracting_GOS()

danmark_strait_coord = np.array([[65,-36],[65,-25]]) #Coord of the West and East limit of danmark strait
Fram_strait_coord = np.array([[79,-17.5],[79,12]]) #Coord of the West and East limit of Fram strait
# Selecting u_gos at Fram Strait
u_gos = u_gos.where((u_gos.longitude < -15) & (u_gos.longitude > -50) & (u_gos.latitude > 60) & (u_gos.latitude < 70))
v_gos = v_gos.where((v_gos.longitude < -15) & (v_gos.longitude > -50) & (v_gos.latitude > 60) & (v_gos.latitude < 70))

current_map(lon, lat, u_gos, v_gos, time,2015,10,15)
