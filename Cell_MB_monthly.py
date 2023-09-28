import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
from datetime import date, timedelta
from scipy import interpolate

#Extracting data
def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc",lat_range = [60,83], lon_range = [-40,20],area_lat = [75,77.5], area_lon = [-10,0]):
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
    
    interp_sit = interpolate.interp1d(time,sit,axis = 0,kind = 'slinear')
    interp_sic = interpolate.interp1d(time,sic,axis = 0,kind = 'slinear')
    sit = interp_sit(np.arange(float(time[0]),float(time[-1])+1))
    sic = interp_sic(np.arange(float(time[0]),float(time[-1])+1))
    time = np.arange(float(time[0]),float(time[-1])+1)
    
    
    #Keeping data only on the area of interest
    
    area_marker = np.where((lat>=area_lat[0]) & (lat<=area_lat[1]) & (lon>=area_lon[0]) & (lon<=area_lon[1]),1,0)
    
    lat = np.where(area_marker == 1,lat,np.nan)
    sit = np.where(area_marker == 1,sit,np.nan)
    sic = np.where(area_marker == 1,sic,np.nan)
    lon = np.where(area_marker == 1,lon,np.nan)
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time, area_marker

def Net_transport_cell(name ='Cell_A',area_lat = [75,77.5], area_lon = [-10,0]):
    """
        This function save as txt in Data/name/ the Mass_bilan of the cell delimited by area_lat and area_lon
    """
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
    month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
    
    lon, lat, sit, sic, sit_uncertainty, time, area_marker = extracting_data_sit(area_lat = area_lat, area_lon = area_lon)
    for year in range(year_,year_end):
        for month in month_names:
            for file in os.listdir(f'Data/X_drift/{year}/{month}'):
                X_drift[f'y{year}'][month].append(np.where(area_marker == 1,np.loadtxt(f'Data/X_drift/{year}/{month}/' + file),np.nan))
                Y_drift[f'y{year}'][month].append(np.where(area_marker == 1,np.loadtxt(f'Data/Y_drift/{year}/{month}/' + file),np.nan))

    northward_transport = np.zeros((year_end -year_,12))
    southward_transport = np.zeros((year_end -year_,12))
    eastward_transport = np.zeros((year_end -year_,12))
    westward_transport = np.zeros((year_end -year_,12))
    Cell_net_transport_monthly = np.zeros((year_end -year_,12))
    for year in range(year_,year_end):
        for month in range(1,13):
            month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
            useful_index = []
            day_of_month_with_sit_data = []
            i = 0
            starting_date = 734419
            for time_ in time:
                corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
                if corresponding_date.year == year and corresponding_date.month == month:
                    useful_index.append(i)
                    day_of_month_with_sit_data.append(corresponding_date.day)
                i += 1
            # All the sit data_array for the month of interest are merged in the following array
            recorded_sit = np.array([sit[n] * sic[n] for n in useful_index])
            recorded_si_drift_X = np.array(X_drift[f"y{year}"][month_names[month-1]][:]) 
            recorded_si_drift_Y = np.array(Y_drift[f"y{year}"][month_names[month-1]][:]) 

            northward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            southward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            eastward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            westward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))

            
            
            for day in range(len(recorded_sit)):
                for line in range(len(recorded_sit[day])):
                    for col in range(len(recorded_sit[day][0])):
                        if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line-1,col]): #in this case we are in the northerest side of the cell
                            northward_si_drift[line,col] = recorded_si_drift_X[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
                        
                        if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line,col-1]): #in this case we are in the westerest side of the cell
                            westward_si_drift[line,col] = recorded_si_drift_X[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360)

                        if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line,col+1]): #in this case we are in the easterest side of the cell
                            eastward_si_drift[line,col] = -(recorded_si_drift_X[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360))
                        
                        if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line+1,col]): #in this case we are in the southerest side of the cell
                            southward_si_drift[line,col] = - recorded_si_drift_X[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) - recorded_si_drift_Y[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
                
                #Transport_ are array covering the same surface as recorded_si_drift and recorded_si_drift. Where there is sea ice drift data, the cell is filled with siv*si_drift.
                northward_transport_current = np.where(abs(northward_si_drift) >0 , northward_si_drift*1000 * recorded_sit[day], np.nan)
                southward_transport_current = np.where(abs(southward_si_drift) >0 , southward_si_drift*1000 * recorded_sit[day], np.nan)
                eastward_transport_current = np.where(abs(eastward_si_drift) >0 , eastward_si_drift*1000 * recorded_sit[day], np.nan)
                westward_transport_current = np.where(abs(westward_si_drift) >0 , westward_si_drift*1000 * recorded_sit[day], np.nan)

                northward_transport[year-year_,month-1] += np.nansum(northward_transport_current * 80000)
                southward_transport[year-year_,month-1] += np.nansum(southward_transport_current * 80000)
                eastward_transport[year-year_,month-1] += np.nansum(eastward_transport_current * 80000)
                westward_transport[year-year_,month-1] += np.nansum(westward_transport_current * 80000)
            #In m^3/month
            Cell_net_transport_monthly[year-year_,month-1] = northward_transport[year-year_,month-1] + southward_transport[year-year_,month-1] + eastward_transport[year-year_,month-1] + westward_transport[year-year_,month-1]
    
    return Cell_net_transport_monthly

def cell_mass_bilan(name ='Cell_A',area_lat = [75,77.5], area_lon = [-10,0]):
    """
        Compute the mass bilan of sea ice melted in the cell during each month by comparing the variation of sea ice along the month with the net export (or import) along the month
    """
    
    lon, lat, sit, sic, sit_uncertainty, time, area_marker = extracting_data_sit(area_lat = area_lat, area_lon = area_lon)
    
    SIV_variation_monthly = np.zeros((year_end -year_,12))
    
    for year in range(year_,year_end):
        for month in range(1,13):
            useful_index = []
            day_of_month_with_sit_data = []
            i = 0
            starting_date = 734419
            for time_ in time:
                corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
                if corresponding_date.year == year and corresponding_date.month == month:
                    useful_index.append(i)
                    day_of_month_with_sit_data.append(corresponding_date.day)
                i += 1
            # All the sit data_array for the month of interest are merged in the following array
            recorded_sit = np.array([sit[n] * sic[n] for n in useful_index])
            SIV_variation_monthly[year-year_,month-1] = (np.nansum(recorded_sit[-1]) -  np.nansum(recorded_sit[0])) * 80000**2 #[m^3], positive when sea ice volume increase over the cell
    Cell_net_transport_monthly = Net_transport_cell(name =name,area_lat = area_lat, area_lon = area_lon) #[m^3] positive when net import of sea ice on the cell
    Sea_ice_melt_on_the_cell = (Cell_net_transport_monthly-SIV_variation_monthly)*1e-9 #[km^3] positive when sea ice melt, negative when sea ice formed
    
    np.savetxt('Data/'+name+'/Cell_'+name[-1]+'_MB_monthly.txt',Sea_ice_melt_on_the_cell)
    return Sea_ice_melt_on_the_cell

year_, year_end = 2011,2020
cell_mass_bilan(name = 'Cell_A',area_lat=[77.5,80],area_lon = [-10,0])
cell_mass_bilan(name = 'Cell_B',area_lat=[75,77.5],area_lon = [-10,0])
cell_mass_bilan(name = 'Cell_C',area_lat=[72.5,75],area_lon = [-20,-10])
cell_mass_bilan(name = 'Cell_D',area_lat=[70,72.5],area_lon = [-20,-10])
""" Sea_ice_melt_on_the_cell = cell_mass_bilan()
print(np.shape(Sea_ice_melt_on_the_cell))
for i in range(len(Sea_ice_melt_on_the_cell)):
    plt.plot([month for month in range(1,13)],Sea_ice_melt_on_the_cell[i])
    plt.grid()
    plt.show() """
