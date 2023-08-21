from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from datetime import datetime, date, timedelta

######### - Function def - ########

def extracting_data(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [64,85], lon_range = [-40,20]):
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
    sic = ds['Sea_Ice_Concentration'].where((ds.Sea_Ice_Thickness != 0) & (ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))    
    sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))
    time =  ds['Time']
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time

def give_time_snapshot(year,month,day):
    """   
        Return the n indix of time[n] given a date. Bug to fix : give_date(give_time_snapshot(date)) doesn't return exactly the date.
    """
    days_from_start = (date(year,month,day) - date(2010,10,1)).days
    time_coverage = time[-1] - time[0] # In days
    time_step = time_coverage/len(time) # In days
    time_snapshot = days_from_start/time_step
    return int(time_snapshot)

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

def mean_SIV(year, month):
    """
        Compute the mean SIV of the month
    """

    # Computing the mensual mean SIT
    useful_date = []
    useful_index = []
    i = 0
    for time_ in time:
        corresponding_date = date(1,1,1) + timedelta(days = int(time_)-365)
        if corresponding_date.year == year and corresponding_date.month == month:
            useful_date.append(int(time_))
            useful_index.append(i)
        i += 1
    recorded_sit = np.array([sit.sel(t = n) * sic.sel(t = n) for n in useful_index])
    mean_SIT = np.nanmean(recorded_sit, axis = 0) 
    # Retrieving the mensual mean SIV
    surface_grid_cell = 80000**2 # [m^2]
    mean_SIV = surface_grid_cell * np.nansum(mean_SIT) # [m^3]

    return mean_SIV
########## - Main - #############

lon, lat, sit, sic, sit_uncertainty, time = extracting_data()


# Each row of SIV_montlhy is a year, each column a month. Their are completed by the total mean Volume value
SIV_monthly = np.zeros((9,12))

for year in range(2011,2020):
    print(f"{year}/2019")
    for month in range(1,13):
        SIV_monthly[year - 2011, month - 1] = mean_SIV(year,month)
month_name = ['january','february','march','april','may','june','july','august','september','october','november','december']
regression_line_slope = []

for month in range(len(SIV_monthly[0])):
    #plt.figure(figsize=(10,5))    

    #liner regression
    b,a = np.polyfit([n for n in range(2011,2020)], [SIV_monthly[y,month] for y in range(len(SIV_monthly))],1)
    regression_line_slope.append(b)
    """plt.plot([n for n in range(2011,2020)],[a+b*n for n in range(2011,2020)],label = "linear regression")
    plt.plot([n for n in range(2011,2020)], [SIV_monthly[y,month] for y in range(len(SIV_monthly))])
    plt.legend()
    plt.grid()
    plt.title(f"Evolution of {month_name[month]} SIV between 2011 and 2019")
    plt.ylabel("SIV [m^3]")
    plt.xlabel("year")
    plt.ylim(0,2*1e12)
    plt.savefig(f"Plots/mensual_evolution/SIV_mensual_{month_name[month]}.png")
    plt.clf()"""

plt.plot([month for month in range(1,13)],[regression_line_slope[month] for month in range(12)])
plt.grid()
plt.title("Slope of the regression line  between 2011 and 2019 for each month")
plt.ylabel("[m^3/year]")
plt.xlabel("month")
plt.savefig(f"Plots/mensual_evolution/SIV_mensual_{month_name[month]}.png")
plt.savefig(f"Plots/mensual_evolution/lin_reg_slope.png")
