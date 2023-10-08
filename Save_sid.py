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
    This script interpol the sea ice drift data over the sit grid and with a daily resolution (the original resolution of this dataset) save it in ./Data
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

    """ plt.subplot(311)
    plt.imshow() """

    #Keeping data only on the gate
    gate_lat = 80
    gate_marker = np.zeros(lat.shape)
    for i in range(len(lat)):
        gate_marker[i,int((abs(lat[i,:]-gate_lat).argmin()))] = 1
    
    """
    lat = lat.where((gate_marker == 1))
    sit = sit.where((gate_marker == 1))
    sic = sic.where((gate_marker == 1))
    lon = lon.where((gate_marker == 1)) """
    gate_marker = 0
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time, gate_marker

def extracting_SI_drift(lat_range = [60,83], lon_range = [-40,20]):
    
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
    for year in range(year_,year_end):
        day_ignored = 0

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
                try:
                    os.makedirs(f"Data/X_drift/{year}/{month_name}")
                except:
                    pass
                try:
                    os.makedirs(f"Data/Y_drift/{year}/{month_name}")
                except:
                    pass
                file =  dir+f"{year}" + "/" + f"{month:02d}" + f"/ice_drift_nh_ease2-750_cdr-v1p0_24h-{year}{month:02d}{day:02d}1200.nc"
                ds = xr.open_dataset(file, decode_times = False)
                
                dX = ds['dX'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                dY = ds['dY'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                
                dX = dX.sel(time = int(ds.time))
                dY = dY.sel(time = int(ds.time))


                if year == year_:
                    if month == 1:
                        if day == 1:
                            lat = ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                            lon = ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)
                            
                ds.close
                
                #Interpolation over SIT grid
                points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
                values_dX = []
                values_dY = []
                for i in range(len(lat)):
                    for j in range(len(lon[0])):
                        if dX[i,j] !=0 and dY[i,j] != 0 and not np.isnan(dX[i,j]) and not np.isnan(dY[i,j]):
                            points.append([lat[i,j],lon[i,j]])
                            values_dX.append(dX[i,j])
                            values_dY.append(dY[i,j])
                points = np.array(points)
                values_dX = np.array(values_dX)
                values_dY = np.array(values_dY)
                if len(values_dX) > 2:
                    dX_interp = interpolate.griddata(points, values_dX, (lat_sit, lon_sit), method='linear')
                    dY_interp = interpolate.griddata(points, values_dY, (lat_sit, lon_sit), method='linear')

                    """ #Keeping data only on the gate
                    dX_interp = np.where(gate_marker == 1, dX_interp,np.nan)
                    dY_interp = np.where(gate_marker == 1, dY_interp,np.nan) """
                    
                    X_drift[f'y{year}'][month_name].append(dX_interp)
                    Y_drift[f'y{year}'][month_name].append(dY_interp)
                else:
                    print('Manual interpolation')
                    #In this case, there are not enough data points (less than 3) to perform a classical interpolation
                    day_ignored += 1
                    manualy_interpolated_X = np.empty(lat_sit.shape)
                    manualy_interpolated_Y = np.empty(lat_sit.shape)
                    manualy_interpolated_X[:] = np.nan
                    manualy_interpolated_Y[:] = np.nan
                    for coords, temp_dX, temp_dY in zip(points,values_dX,values_dY):
                        idx_lat = int((np.abs(lat_sit - coords[0])).argmin())
                        idx_lon = int((np.abs(lon_sit - coords[1])).argmin())
                        
                        
                        manualy_interpolated_X[np.unravel_index(idx_lon,lon_sit.shape)] = temp_dX
                        manualy_interpolated_Y[np.unravel_index(idx_lon,lon_sit.shape)] = temp_dY

                    """ #Keeping data only on the gate
                    manualy_interpolated_X = np.where(gate_marker == 1, manualy_interpolated_X,np.nan)
                    manualy_interpolated_Y = np.where(gate_marker == 1, manualy_interpolated_Y,np.nan) """

                    X_drift[f'y{year}'][month_name].append(manualy_interpolated_X)
                    Y_drift[f'y{year}'][month_name].append(manualy_interpolated_Y)
                try:
                    os.remove(f'Data/X_drift/{year}/{month_name}/{day}.txt')
                except:
                    pass
                try:
                    os.remove(f'Data/Y_drift/{year}/{month_name}/{day}.txt')
                except:
                    pass
                np.savetxt(f'Data/X_drift/{year}/{month_name}/{day:02d}.txt',X_drift[f'y{year}'][month_name][-1])
                np.savetxt(f'Data/Y_drift/{year}/{month_name}/{day:02d}.txt',Y_drift[f'y{year}'][month_name][-1])
        print(f"{day_ignored} days man interp during year {year}")
    
    return X_drift, Y_drift, lat, lon



if __name__ == '__main__':
    
    year_ = 2011
    year_end = 2021
    #X_drift, Y_drift, lat, lon = extracting_SI_drift()
    lon_sit,lat_sit, sit, sic, sit_uncertainty, sit_time, gate_marker= extracting_data_sit()
    X_drift, Y_drift, lat_drift, lon_drift = extracting_SI_drift()
    