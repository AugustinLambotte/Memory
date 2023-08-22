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


def extracting_SI_drift(lat_range = [79,80], lon_range = [-40,20]):
    """
        Extract the SI drift data at Fram Strait
    """
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
                dX = ds['dX'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat  > lat_min) & (ds.lat < lat_max), drop = True)
                dY = ds['dY'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat  > lat_min) & (ds.lat < lat_max), drop = True)
                
                dX = dX.sel(time = int(ds.time))
                dY = dY.sel(time = int(ds.time))

                #dX = dX.dropna(dim = "yc",how = "all")
                #dY = dY.dropna(dim = "yc",how = "all")

                if year == 2010:
                    if month == 1:
                        if day == 1:
                            lat = ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat  > lat_min) & (ds.lat < lat_max), drop = True)
                            lon = ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat  > lat_min) & (ds.lat < lat_max), drop = True)
                            xc = ds['xc']
                            yc = ds['yc']
                ds.close
                
                X_drift[f'y{year}'][month_name].append(dX)
                Y_drift[f'y{year}'][month_name].append(dY)
    
    
    return X_drift, Y_drift, lat, lon, xc, yc

def extracting_sit(lat_range = [79,80], lon_range = [-40,20]):
    """
        Extract the sit data at Fram Strait
    """
    file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc"

    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_min = lat_range[0]
    lat_max = lat_range[1]

    ds = xr.open_dataset(file, decode_times = False)
    lon = ds['Longitude'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude < lat_max) & (ds.Latitude > lat_min), drop = True)
    lat = ds['Latitude'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude < lat_max) & (ds.Latitude > lat_min), drop = True)
    sit = ds['Sea_Ice_Thickness'].where((ds.Sea_Ice_Thickness != 0) & (ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude < lat_max) & (ds.Latitude > lat_min), drop = True)
    sic = ds['Sea_Ice_Concentration'].where((ds.Sea_Ice_Thickness != 0) & (ds.Longitude > lon_min) & (ds.Longitude < lon_max) &(ds.Latitude < lat_max) & (ds.Latitude > lat_min), drop = True)
    time = ds['Time']
    sit = sit*sic
    ds.close
    return lon, lat, sit, time

def Fram_strait_mass_bilan(year, month):
    """
        Return the monthly mean SIV that pass through the Fram strait during the given month
    """
    print(year, month)
    month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
    useful_date = []
    useful_index = []
    day_of_month_with_sit_data = []
    i = 0
    starting_date = 734419
    for time_ in sit_time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
        if corresponding_date.year == year and corresponding_date.month == month:
            useful_date.append(int(time_))
            useful_index.append(i)
            day_of_month_with_sit_data.append(corresponding_date.day)
        i += 1
    # All the sit data_array for the month of interest are merged in the following array
    recorded_sit = np.array([sit.sel(t = n) for n in useful_index])
    recorded_si_drift = np.array([X_drift[f"y{year}"][month_names[month-1]][day-1] for day in day_of_month_with_sit_data])
    if len(recorded_si_drift) == 0 or len(recorded_sit) == 0:
        return 0
    
    # The following lines construct lat_with_drift and lon_with drift. Two array with dimension 
    # [day_in_the_month_with_sit_data, nb_cell_with_drift_data]. Each element are the latitude and longitude of a cell with sea ice drift
    # during a day when we have a record of sit.
    drift_col = []
    drift_line = []
    lat_with_drift = []
    lon_with_drift = []
    SI_drift = []
    for day in range(len(recorded_si_drift)):
        drift_col.append([])
        drift_line.append([])
        lat_with_drift.append([])
        lon_with_drift.append([])
        SI_drift.append([])
        for line in range(len(recorded_si_drift[day])):
            for col in range(len(recorded_si_drift[day,line,:])):
                if not np.isnan(recorded_si_drift[day,line,col]):
                    drift_line[-1].append(line)
                    drift_col[-1].append(col)
        for line,col in zip(drift_line[day],drift_col[day]):
            lat_with_drift[-1].append(float(drift_lat[line,col]))
            lon_with_drift[-1].append(float(drift_lon[line,col]))
            SI_drift[-1].append(float(recorded_si_drift[day,line,col]))

    
    # Same as above but for lat and lon with sit data
    sit_col = []
    sit_line = []
    lat_with_sit = []
    lon_with_sit = []
    SIT = []
    for day in range(len(recorded_sit)):
        sit_col.append([])
        sit_line.append([])
        lat_with_sit.append([])
        lon_with_sit.append([])
        SIT.append([])
        for line in range(len(recorded_sit[day])):
            for col in range(len(recorded_sit[day,line,:])):
                if not np.isnan(recorded_sit[day,line,col]):
                    sit_line[-1].append(line)
                    sit_col[-1].append(col)
        for line,col in zip(sit_line[day],sit_col[day]):
            lat_with_sit[-1].append(float(sit_lat[line,col]))
            lon_with_sit[-1].append(float(sit_lon[line,col]))
            SIT[-1].append(float(recorded_sit[day,line,col]))
            
      
    # We now pass over the longitude where we have sit data and we search the closer (in longitude) drift data.
    # We stock in sit_drift tulpes with the cell SIT and the corresponding (i.e. closer) drift value

    sit_drift = []

    for day in range(len(recorded_si_drift)):
        sit_drift.append([])
        for longitude_sit in lon_with_sit[day]:
            index = np.argmin([abs(longitude_drift - longitude_sit) for longitude_drift in lon_with_drift[day]])
            sit_drift[-1].append(dict(sit = SIT[day][index], drift = SI_drift[day][index]))
    
    #Now we compute the volume of ice passing each of the day where we have data. Note that the SIT data are upon a grid_cell of 80km^2
    SIV_passing = np.zeros(len(recorded_si_drift))
    for day in range(len(recorded_si_drift)):
        daily_volume = 0
        for cell in (sit_drift[day]):
            daily_volume += cell['sit'] * cell['drift']*1000 * 80000 #in m^3
        SIV_passing[day] = daily_volume

    monthly_averaged_SIV = np.mean(SIV_passing)
    
    return monthly_averaged_SIV

##### - Main - ####

X_drift, Y_drift, drift_lat, drift_lon, xc, yc = extracting_SI_drift()
sit_lon, sit_lat, sit, sit_time = extracting_sit()

for year in range(2010,2021):
    Mass_bilan_fram_strait = np.zeros(12)
    for month in range(1,13):
        Mass_bilan_fram_strait[month-1] = Fram_strait_mass_bilan(year,month)

    plt.plot([month for month in range(1,13)], Mass_bilan_fram_strait)
    plt.title(f'monthly averaged daily SIV passing through Fram Strait Northward in [m^3] - {year}')
    plt.grid()
    plt.savefig(f"Plots/Fram_strait/{year}-Mass_bilan_Fram_strait.png")
    plt.clf()

