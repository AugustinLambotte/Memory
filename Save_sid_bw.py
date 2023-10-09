import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime, date, timedelta
import matplotlib.path as mpath
import cmocean as cmo
import cartopy.crs as ccrs
import scipy
import os 
"""
    This script interpol the sea ice drift data over the sit grid with a bi-weekly resolution (the temporal resolution of cryosat) and save it in ./Data
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


if __name__ == '__main__':
    
    year_ = 2010
    year_end = 2021
    lon_sit,lat_sit, sit, sic, sit_uncertainty, sit_time, gate_marker= extracting_data_sit()

    #Day flatten sit is used in the temporal interpolation of sid over the bi-weekly grid.
    day_flatten_sit = []
    starting_date = 734419
    #first_jan_day = starting_date + (date(2011,1,1) - date(2010,10,1)).days #The 1st january 2011 in the original sit dataset date format
    rec_date = []
    for time_ in sit_time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
        if corresponding_date.year >= year_ and corresponding_date.year < year_end:
            day_flatten_sit.append(int(time_)-starting_date)
            rec_date.append(corresponding_date)
            try:
                    os.makedirs(f"Data/bw/X_drift")
                    os.makedirs(f"Data/bw/Y_drift")
            except:
                pass
    
    X_drift_flatten = []
    Y_drift_flatten = []
    day_flatten_sid = []
    day = 0
    for year in range(year_,year_end):
        if year == 2010:
            month_names = ["oct","nov", "dec"]
        else:
            month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]

        for month in month_names:
            for file in os.listdir(f'Data/X_drift/{year}/{month}'):
                day_flatten_sid.append(day)
                X_drift_flatten.append(np.loadtxt(f'Data/X_drift/{year}/{month}/' + file))
                Y_drift_flatten.append(np.loadtxt(f'Data/Y_drift/{year}/{month}/' + file))
                day += 1
    #Temporal interpolation. We interpolate this daily series over the biweekly SIV time resolution
    x_interpolator = interpolate.interp1d(day_flatten_sid,X_drift_flatten,axis = 0)
    y_interpolator = interpolate.interp1d(day_flatten_sid,Y_drift_flatten,axis = 0)

    X_interpolated = x_interpolator(day_flatten_sit)
    Y_interpolated = y_interpolator(day_flatten_sit)
    for date_, X, Y in zip(rec_date,X_interpolated,Y_interpolated):
        np.savetxt(f'Data/bw/X_drift/{date_}.txt',X)
        np.savetxt(f'Data/bw/Y_drift/{date_}.txt',Y)
    #The following array save the date where we have data in two format. yyyy-mm-dd and day from 2011-01-01

    saving_date = np.zeros((len(day_flatten_sit),4))
    for i in range(len(day_flatten_sit)):
        saving_date[i,0] = rec_date[i].year
        saving_date[i,1] = rec_date[i].month
        saving_date[i,2] = rec_date[i].day
        saving_date[i,3] = day_flatten_sit[i]
    np.savetxt('Data/bw/date.txt',saving_date, fmt='%i')
    #X_drift, Y_drift, lat_drift, lon_drift = extracting_SI_drift()
    