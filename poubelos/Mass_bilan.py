from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from datetime import datetime, date, timedelta

from Useful_function import extracting_SI_velocity, give_date, give_time_snapshot


def SI_vel_map(lon, lat, usi, vsi, year, month, day, projection = ccrs.LambertConformal(central_longitude = -20), figsize = (14,10)):
    """
        Plot Sea ice velocity given lon,lat, usi, vsi, year, month and day
    """
    time_snapshot = (date(year,month,day) - date(2010,10,1)).days
    current_magnitude = np.sqrt(usi**2 + vsi**2)
    fig,axs = plt.subplots(figsize = figsize, subplot_kw={'projection': projection})
    axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
    axs.coastlines()

    #Magnitude plot
    cs = axs.contourf(lon,lat,current_magnitude[time_snapshot],transform =ccrs.PlateCarree())
    
    #Vector plot
    mymap = plt.streamplot(np.array(lon[:]),np.array(lat[:]),np.array(usi[time_snapshot]),np.array(vsi[time_snapshot]), transform=ccrs.PlateCarree(), density=4)
    #axs.quiver(np.array(lon[:]),np.array(lat[:]),np.array(u_gos[time_snapshot]),np.array(v_gos[time_snapshot]), transform = ccrs.PlateCarree())
    fig.colorbar(cs, ax=axs)
    plt.title(f"{year} - {month} - {day} : Sea Ice velocity [m/s]")
    plt.show()

def mass_bilan_Fram_Strait(vsi, year, month, day):
    #usi = usi.where((usi.longitude < -15) & (usi.longitude > -50) & (usi.latitude == 78))
    vsi = vsi.where((vsi.longitude < 12) & (vsi.longitude > -16) & (vsi.latitude == 78))

    time_snapshot = (date(year,month,day) - date(2010,10,1)).days
    velocity = []
    for day in range(len(vsi.time)):
        velocity.append(np.nanmean(vsi.isel(time = day)))
    plt.plot([day for day in range(len(vsi.time))], velocity)
    plt.title("Sea ice velocity at fram strait day to day [m/s]")
    plt.grid()
    plt.xlabel("days from 2010-10-01")
    plt.ylabel("Sea ice velocity [m/s]")
    plt.show()

def mass_bilan_Denmark_Strait(vsi):
    #usi = usi.where((usi.longitude < -15) & (usi.longitude > -50) & (usi.latitude == 78))
    vsi = vsi.where((vsi.longitude < -22) & (vsi.longitude > -34.6) & (vsi.latitude == 66.5))

    #time_snapshot = (date(year,month,day) - date(2010,10,1)).days
    velocity = []
    for day in range(len(vsi.time)):
        velocity.append(np.nanmean(vsi.isel(time = day)))
    plt.plot([day for day in range(len(vsi.time))], velocity)
    plt.title("Sea ice velocity at Denamrk strait day to day [m/s]")
    plt.grid()
    plt.xlabel("days from 2010-10-01")
    plt.ylabel("Sea ice velocity [m/s]")
    plt.show()

lon,lat,usi,vsi, time = extracting_SI_velocity()
#mass_bilan_Fram_Strait(vsi,2012,9,15)
#SI_vel_map(lon,lat,usi,vsi,2012,3,15)
mass_bilan_Denmark_Strait(vsi)

