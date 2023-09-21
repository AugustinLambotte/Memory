from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import xarray as xr
import pandas as pd
from scipy import interpolate
from datetime import datetime, date, timedelta
import cmocean
from SIT_map import extracting_data_sit
year_ = 2010
year_end = 2021
def extracting_SI_drift(lat_range = [64,80.5], lon_range = [-40,20]):
    
    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    #print(ds.Latitude.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max), drop = True))

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
    for year in range(year_,year_end):
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
                
                dX = ds['dX'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max) & (ds.lat > 65.4 + (76.5-65.4)/(9+17) * (ds.lon + 17)), drop = True)
                dY = ds['dY'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max) & (ds.lat > 65.4 + (76.5-65.4)/(9+17) * (ds.lon + 17)), drop = True)
                
                dX = dX.sel(time = int(ds.time))
                dY = dY.sel(time = int(ds.time))


                if year == year_:
                    if month == 1:
                        if day == 1:
                            lat = ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max) & (ds.lat > 65.4 + (76.5-65.4)/(9+17) * (ds.lon + 17)), drop = True)
                            lon = ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max) & (ds.lat > 65.4 + (76.5-65.4)/(9+17) * (ds.lon + 17)), drop = True)
                ds.close
                
                X_drift[f'y{year}'][month_name].append(dX)
                Y_drift[f'y{year}'][month_name].append(dY)
    
    return X_drift, Y_drift, lat, lon

def mensual_mean(year, month, interpolated = False):
    month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
    print(year,month)
    X_drift_averaged = np.nanmean(X_drift[f"y{year}"][f"{month_names[month-1]}"], axis = 0)
    Y_drift_averaged = np.nanmean(Y_drift[f"y{year}"][f"{month_names[month-1]}"], axis = 0)

    
    if interpolated:
        #Interpolation over SIT grid
        points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
        values_dX = []
        values_dY = []
        for i in range(len(lat_drift)):
            for j in range(len(lon_drift[0])):
                if X_drift_averaged[i,j] !=0 and Y_drift_averaged[i,j] != 0 and not np.isnan(X_drift_averaged[i,j]) and not np.isnan(Y_drift_averaged[i,j]):
                    points.append([lat_drift[i,j],lon_drift[i,j]])
                    values_dX.append(X_drift_averaged[i,j])
                    values_dY.append(Y_drift_averaged[i,j])
        points = np.array(points)
        values_dX = np.array(values_dX)
        values_dY = np.array(values_dY)
        dX_interp = interpolate.griddata(points, values_dX, (lat_sit, lon_sit), method='linear')
        dY_interp = interpolate.griddata(points, values_dY, (lat_sit, lon_sit), method='linear')
        
        return dX_interp, dY_interp, lon_sit, lat_sit
    else:
        return X_drift_averaged, Y_drift_averaged, lon_drift, lat_drift

if __name__ == "__main__":

    X_drift, Y_drift, lat_drift, lon_drift = extracting_SI_drift()
    lon_sit,lat_sit,_,_,_,_ = extracting_data_sit()
    interpolated = True
    for year in range(year_,year_end):
        for month in range(1,13):
            X_drift_av, Y_drift_av, lon,lat = mensual_mean(year,month,interpolated=interpolated)
            current_magnitude = np.sqrt(X_drift_av**2 + Y_drift_av**2)
            current_magnitude_ms = np.sqrt(X_drift_av**2 + Y_drift_av**2) * 1000/(24*60*60) #Passing from km/day in m/s
            
            #fig = plt.figure(figsize=(10,10))
            #axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
            
            fig,axs = plt.subplots(figsize = (10,10), subplot_kw={'projection': ccrs.LambertConformal(central_longitude = -18)})
            #axs.set_extent([-40, 11, 62, 85], crs = ccrs.PlateCarree())
            
            xlim = [-43, 16]
            ylim = [61, 81]
            lower_space = 3 
            rect = mpath.Path([[xlim[0], ylim[0]],
                            [xlim[1], ylim[0]],
                            [xlim[1], ylim[1]],
                            [xlim[0], ylim[1]],
                            [xlim[0], ylim[0]],
                            ]).interpolated(20)
            proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
            rect_in_target = proj_to_data.transform_path(rect)
            axs.set_boundary(rect_in_target)
            axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
            axs.coastlines()
            axs.gridlines()
            axs.set_title(f"Sea Ice mean drift in [m/s] during {month}/{year}")
            #Magnitude plot
            levels = np.linspace(0,0.3,10)
            cs = axs.contourf(lon,lat,current_magnitude_ms, cmap = "cmo.speed", levels = levels, transform =ccrs.PlateCarree())
            #Vector plot
            X_drift_av = X_drift_av/current_magnitude * 10
            Y_drift_av = Y_drift_av/current_magnitude * 10
            axs.quiver(np.array(lon),np.array(lat),np.array(X_drift_av),np.array(Y_drift_av),scale = 400,transform = ccrs.PlateCarree())
            fig.colorbar(cs, ax=axs,ticks = [0,0.1,0.2,0.3])
            if interpolated:
                plt.savefig(f"Plots/mean/Sea_ice_drift/interpolated/{year}/SI_drift_{year}-{month}.png")
            else:
                plt.savefig(f"Plots/mean/Sea_ice_drift/Original/{year}/SI_drift_{year}-{month}.png")

            plt.clf()
    
                
                

