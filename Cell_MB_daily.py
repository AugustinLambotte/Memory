import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
from datetime import date, timedelta
from scipy import interpolate,fft


#Extracting data
def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", daily_interpolation = False,lat_range = [60,83], lon_range = [-40,20],area_lat = [75,77.5], area_lon = [-10,0]):
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

    #Interpolation to daily resolution
    interp_sit = interpolate.interp1d(time,sit,axis = 0,kind = 'linear', assume_sorted = True)
    interp_sic = interpolate.interp1d(time,sic,axis = 0,kind = 'linear', assume_sorted = True)
    sit = interp_sit(np.arange(float(time[0]),float(time[-1])+1))
    sic = interp_sic(np.arange(float(time[0]),float(time[-1])+1))
    time = np.arange(float(time[0]),float(time[-1])+1)
    

    #Keeping data only on the area of interest
    
    area_marker = np.where((lat>=area_lat[0]) & (lat<=area_lat[1]) & (lon>=area_lon[0]) & (lon<=area_lon[1]),1,np.nan)
    
    """ lat = np.where(area_marker == 1,lat,np.nan)
    sit = np.where(area_marker == 1,sit,np.nan)
    sic = np.where(area_marker == 1,sic,np.nan)
    lon = np.where(area_marker == 1,lon,np.nan) """

    lat = np.array(lat)
    sit = np.array(sit)
    sic = np.array(sic)
    lon = np.array(lon)
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
                X_drift[f'y{year}'][month].append(np.loadtxt(f'Data/X_drift/{year}/{month}/' + file))
                Y_drift[f'y{year}'][month].append(np.loadtxt(f'Data/Y_drift/{year}/{month}/' + file))
    
    Cell_net_transport_daily = []
    for year in range(year_,year_end):
        Cell_net_transport_daily.append([])
        useful_index = []
        day_of_year_with_sit_data = []
        i = 0
        starting_date = 734419
        for time_ in time:
            corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
            if corresponding_date.year == year:
                
                useful_index.append(i)
                day_of_year_with_sit_data.append(corresponding_date.day)
            i += 1
        # All the sit data_array for the month of interest are merged in the following array
        recorded_sit = np.array([sit[n] * sic[n] for n in useful_index])
        recorded_si_drift_X = np.zeros(np.shape(recorded_sit))
        recorded_si_drift_Y = np.zeros(np.shape(recorded_sit))
        j=0
        for month in range(12):
            i = 0
            while i != len(X_drift[f"y{year}"][month_names[month]][:]):
                #print(X_drift[f"y{year}"][month_names[month]][i])
                recorded_si_drift_X[j] = X_drift[f"y{year}"][month_names[month]][i] 
                recorded_si_drift_Y[j] = Y_drift[f"y{year}"][month_names[month]][i] 
                i += 1
                j += 1
        for day in range(len(recorded_sit)):
            northward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            southward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            eastward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            westward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
            
            
            for line in range(len(recorded_sit[day])-1):
                for col in range(len(recorded_sit[day][0])-1): #The minus one for col and line are to handle "effet de bords" when accessing area_marker[line+1,col] for example: ATTENTION it will cause lack of data if area marker is on the edge of the map
                    if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line+1,col]): #in this case we are in the northerest side of the cell
                        northward_si_drift[line,col] = recorded_si_drift_X[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) - recorded_si_drift_Y[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
                    
                    if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line,col+1]): #in this case we are in the westerest side of the cell
                        westward_si_drift[line,col] = recorded_si_drift_X[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360)

                    if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line,col-1]): #in this case we are in the easterest side of the cell
                        eastward_si_drift[line,col] = -(recorded_si_drift_X[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360))
                    
                    if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line-1,col]): #in this case we are in the southerest side of the cell
                        southward_si_drift[line,col] = - recorded_si_drift_X[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
            #Transport_ are array covering the same surface as recorded_si_drift and recorded_si_drift. Where there is sea ice drift data, the cell is filled with siv*si_drift in meters.
            northward_transport_current = np.where(abs(northward_si_drift) >0 , northward_si_drift*1000 * recorded_sit[day], np.nan)
            southward_transport_current = np.where(abs(southward_si_drift) >0 , southward_si_drift*1000 * recorded_sit[day], np.nan)
            eastward_transport_current = np.where(abs(eastward_si_drift) >0 , eastward_si_drift*1000 * recorded_sit[day], np.nan)
            westward_transport_current = np.where(abs(westward_si_drift) >0 , westward_si_drift*1000 * recorded_sit[day], np.nan)
            
            
            Cell_net_transport_daily[-1].append(np.nansum(northward_transport_current * 80000)+ np.nansum(southward_transport_current * 80000) + np.nansum(eastward_transport_current * 80000) + np.nansum(westward_transport_current * 80000))
            
    return Cell_net_transport_daily

def cell_mass_bilan(name ='Cell_A',area_lat = [75,77.5], area_lon = [-10,0]):
    """
        Compute the mass bilan of sea ice melted in the cell during each month by comparing the variation of sea ice along the month with the net export (or import) along the month
    """
    
    lon, lat, sit, sic, sit_uncertainty, time, area_marker = extracting_data_sit(area_lat = area_lat, area_lon = area_lon)
    
    SIV_variation_daily = []
    cell_net_transport_daily = Net_transport_cell(name =name,area_lat = area_lat, area_lon = area_lon) #[m^3] positive when net import of sea ice on the cell
    Sea_ice_melt_on_the_cell = []
    SIV_daily = []
    for year in range(year_,year_end):
        SIV_daily.append([])
        Sea_ice_melt_on_the_cell.append([])
        SIV_variation_daily.append([])
        useful_index = []
        day_of_year_with_sit_data = []
        i = 0
        starting_date = 734419
        for time_ in time:
            corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
            if corresponding_date.year == year:
                if len(useful_index) == 0:
                    #Before taking the first day of the year we take the last day of the previous year in order to take the net SIV change for the first day of the year.
                    useful_index.append(i-1)
                useful_index.append(i)
                day_of_year_with_sit_data.append(corresponding_date.day)
            i += 1
        # All the sit data_array for the month of interest are merged in the following array
        recorded_sit = np.array([np.where(area_marker == 1, sit[n] * sic[n],np.nan) for n in useful_index])
        for day in range(1,len(recorded_sit)):
            SIV_variation = (np.nansum(recorded_sit[day])-np.nansum(recorded_sit[day-1])) * 80000**2 #[m^3] positive when sea ice volume increase over the cell
            
            Sea_ice_melt_on_the_cell[-1].append((cell_net_transport_daily[year-year_][day-1] - SIV_variation)*1e-9) #[km^3] positive when sea ice melt, negative when sea ice form (could be seen as fresh water flux through water)
            #print(Sea_ice_melt_on_the_cell[-1][-1])
            SIV_variation_daily[-1].append(SIV_variation * 1e-9)#[km^3]
            SIV_daily[-1].append(np.nansum(recorded_sit[day]*80000**2 * 1e-9)) #[km^3]
    
    """ plt.subplot(311)
    plt.plot(range(len(SIV_daily[-1])),SIV_daily[-1])
    plt.title('siv')
    plt.grid()
    plt.subplot(312)
    plt.plot(range(len(Sea_ice_melt_on_the_cell[-1])),Sea_ice_melt_on_the_cell[-1])
    plt.grid()
    plt.title('melting')
    plt.subplot(313)
    plt.grid()
    plt.plot(range(len(cell_net_transport_daily[-1])),cell_net_transport_daily[-1])
    plt.title('transport')
    plt.show() """
    #Turn in format with regular shape numpy array in order to be able to save it
    saving_format_melt = np.zeros((year_end-year_,366))
    saving_format_melt[:] = np.nan
    saving_format_siv = np.zeros((year_end-year_,366))
    saving_format_siv[:] = np.nan
    saving_format_transp = np.zeros((year_end-year_,366))
    saving_format_transp[:] = np.nan
    for year in range(year_,year_end):
        for day in range(len(Sea_ice_melt_on_the_cell[year-year_])):
            saving_format_melt[year-year_,day] = Sea_ice_melt_on_the_cell[year-year_][day]
            saving_format_siv[year-year_,day] = SIV_daily[year-year_][day]
            saving_format_transp[year-year_,day] = cell_net_transport_daily[year-year_][day] *1e-9

            
    np.savetxt('Data/'+name+'/Cell_'+name[-1]+'_fw_daily.txt',saving_format_melt)
    np.savetxt('Data/'+name+'/Cell_'+name[-1]+'_siv_daily.txt',saving_format_siv)
    np.savetxt('Data/'+name+'/Cell_'+name[-1]+'_transp_daily.txt',saving_format_transp)
    return Sea_ice_melt_on_the_cell
year_, year_end = 2011,2020
cell_mass_bilan(name = 'Cell_A',area_lat=[77.5,80],area_lon = [-10,0])
cell_mass_bilan(name = 'Cell_B',area_lat=[75,77.5],area_lon = [-10,0])
cell_mass_bilan(name = 'Cell_C',area_lat=[72.5,75],area_lon = [-20,-10])
cell_mass_bilan(name = 'Cell_D',area_lat=[70,72.5],area_lon = [-20,-10])

""" Sea_ice_melt_on_the_cell = cell_mass_bilan()
print(np.shape(Sea_ice_melt_on_the_cell))
plt.plot(np.arange(len(Sea_ice_melt_on_the_cell[0])),Sea_ice_melt_on_the_cell[0])
plt.grid()
plt.show() """
