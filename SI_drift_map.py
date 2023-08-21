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


def extracting_SI_drift(lat_range = [64,85], lon_range = [-40,20]):
    
    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]
    X_drift = dict( y2010 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2011 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2012 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2013 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2014 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2015 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2016 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2017 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2018 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2019 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2020 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),)
    Y_drift = dict( y2010 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2011 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2012 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2013 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2014 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2015 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2016 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2017 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2018 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2019 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),
                    y2020 = dict(jan = [], feb = [], mar = [], april = [], may = [], june = [], july = [], aug = [], sept = [], oct = [], nov =[], dec = []),)
    dir = "C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/"
    print("\n##### - Extracting Sea Ice drift data - #####\n")
    for year in range(2010,2021):
        for month in range(1,13):
            print(f'{month}/{year}')
            if month ==1:
                nb_days = 31
                month_name = "jan"
            if month ==2:
                month_name = "feb"
                if year == 2010:
                    nb_days = 28
                if year == 2011:
                    nb_days = 28
                if year == 2012:
                    nb_days = 29
                if year == 2013:
                    nb_days = 28
                if year == 2014:
                    nb_days = 28
                if year == 2015:
                    nb_days = 28
                if year == 2016:
                    nb_days = 29
                if year == 2017:
                    nb_days = 28
                if year == 2018:
                    nb_days = 28
                if year == 2019:
                    nb_days = 28
                if year == 2020:
                    nb_days = 29
            if month ==3:
                nb_days = 31
                month_name = "mar"
            if month ==4:
                nb_days = 30
                month_name = "april"
            if month ==5:
                nb_days = 31
                month_name = "may"
            if month ==6:
                nb_days = 30
                month_name = "june"
            if month ==7:
                nb_days = 31
                month_name = "july"
            if month ==8:
                nb_days = 31
                month_name = "aug"
            if month ==9:
                nb_days = 30
                month_name = "sept"
            if month ==10:
                nb_days = 31
                month_name = "oct"
            if month ==11:
                nb_days = 30
                month_name = "nov"
            if month ==12:
                nb_days = 31
                month_name = "dec"
            for day in range(1,nb_days+1):
                file =  dir+f"{year}" + "/" + f"{month:02d}" + f"/ice_drift_nh_ease2-750_cdr-v1p0_24h-{year}{month:02d}{day:02d}1200.nc"
                ds = xr.open_dataset(file, decode_times = False)
                dX = ds['dX'].where((ds.lon1 > lon_min) & (ds.lon1 < lon_max) & (ds.lat1 > lat_min) & (ds.lat1 < lat_max))
                dY = ds['dY'].where((ds.lon1 > lon_min) & (ds.lon1 < lon_max) & (ds.lat1 > lat_min) & (ds.lat1 < lat_max))
                
                dX = dX.sel(time = int(ds.time))
                dY = dY.sel(time = int(ds.time))
                if year == 2010:
                    if month == 1:
                        if day == 1:
                            lat = ds['lat']
                            lon = ds['lon']
                ds.close
                X_drift[f'y{year}'][month_name].append(dX)
                Y_drift[f'y{year}'][month_name].append(dY)
    return X_drift, Y_drift, lat, lon

def mensual_mean(year, month):
    month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
    X_drift_averaged = np.nanmean(X_drift[f"y{year}"][f"{month_names[month-1]}"], axis = 0)
    Y_drift_averaged = np.nanmean(Y_drift[f"y{year}"][f"{month_names[month-1]}"], axis = 0)
    
    return X_drift_averaged, Y_drift_averaged

##### - Main - ####

X_drift, Y_drift, lat, lon = extracting_SI_drift()

for year in range(2010,2021):
    for month in range(1,13):
        X_drift_av, Y_drift_av = mensual_mean(year,month)
        current_magnitude = np.sqrt(X_drift_av**2 + Y_drift_av**2)
        fig,axs = plt.subplots(figsize = (10,10), subplot_kw={'projection': ccrs.LambertConformal(central_longitude = -20)})
        axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
        axs.coastlines()
        axs.set_title(f"Sea Ice mean drift in km/days during {month}/{year}")
        #Magnitude plot
        levels = np.linspace(0,35,30)
        cs = axs.contourf(lon,lat,current_magnitude, cmap = "cmo.speed",levels = levels, transform =ccrs.PlateCarree())
            
        #Vector plot
        axs.quiver(np.array(lon),np.array(lat),np.array(X_drift_av),np.array(Y_drift_av),scale = 400,transform = ccrs.PlateCarree())
        fig.colorbar(cs, ax=axs)
        plt.savefig(f"Plots/mean/Sea_ice_drift/{year}/SI_drift_{year}-{month}.png")
        plt.clf()

             
            
