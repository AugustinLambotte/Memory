import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
from datetime import date, timedelta
from scipy import interpolate,fft

"""
    This code is the same as cell_MB_daily but instead of using interpolated sit over the daily sid data we use interpolated sid data over the sit grid (which is bi-weekely (bw)).
"""
#Extracting data
def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [60,83], lon_range = [-40,20],area_lat = [75,77.5], area_lon = [-10,0]):
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
    
    month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
    X_drift = []
    Y_drift = []
    sid_date = []
    lon, lat, sit, sic, sit_uncertainty, time, area_marker = extracting_data_sit(area_lat = area_lat, area_lon = area_lon)
    for file in os.listdir(f'Data/bw/X_drift'):
        X_drift.append(np.loadtxt(f'Data/bw/X_drift/{file}'))
        Y_drift.append(np.loadtxt(f'Data/bw/Y_drift/{file}'))
        sid_date.append(file)
    
    Cell_net_transport_bw = []
    useful_index = []
    day_of_year_with_sit_data = []
    i = 0
    starting_date = 734419

    for time_ in time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
        if corresponding_date.year >= year_ and corresponding_date.year < year_end:
            useful_index.append(i)
            day_of_year_with_sit_data.append(corresponding_date.day)
        i += 1

    # All the sit data_array for the month of interest are merged in the following array
    recorded_sit = np.array([sit[n] * sic[n] for n in useful_index])
    for day in range(len(recorded_sit)):
        northward_si_drift = np.zeros(np.shape(X_drift[0]))
        southward_si_drift = np.zeros(np.shape(X_drift[0]))
        eastward_si_drift = np.zeros(np.shape(X_drift[0]))
        westward_si_drift = np.zeros(np.shape(X_drift[0]))
        
        
        for line in range(len(recorded_sit[day])-1):
            for col in range(len(recorded_sit[day][0])-1): #The minus one for col and line are to handle "effet de bords" when accessing area_marker[line+1,col] for example: ATTENTION it will cause lack of data if area marker is on the edge of the map
                if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line+1,col]): #in this case we are in the northerest side of the cell
                    northward_si_drift[line,col] = X_drift[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) - Y_drift[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
                
                if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line,col+1]): #in this case we are in the westerest side of the cell
                    westward_si_drift[line,col] = X_drift[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + Y_drift[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360)

                if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line,col-1]): #in this case we are in the easterest side of the cell
                    eastward_si_drift[line,col] = -(X_drift[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + Y_drift[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360))
                
                if np.isnan(area_marker[line,col]) and not np.isnan(area_marker[line-1,col]): #in this case we are in the southerest side of the cell
                    southward_si_drift[line,col] = - X_drift[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) + Y_drift[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
        #Transport_ are array covering the same surface as recorded_si_drift and recorded_si_drift. Where there is sea ice drift data, the cell is filled with siv*si_drift in meters.
        northward_transport_current = np.where(abs(northward_si_drift) >0 , northward_si_drift*1000 * recorded_sit[day], np.nan)
        southward_transport_current = np.where(abs(southward_si_drift) >0 , southward_si_drift*1000 * recorded_sit[day], np.nan)
        eastward_transport_current = np.where(abs(eastward_si_drift) >0 , eastward_si_drift*1000 * recorded_sit[day], np.nan)
        westward_transport_current = np.where(abs(westward_si_drift) >0 , westward_si_drift*1000 * recorded_sit[day], np.nan)
        
        
        Cell_net_transport_bw.append(np.nansum(northward_transport_current * 80000)+ np.nansum(southward_transport_current * 80000) + np.nansum(eastward_transport_current * 80000) + np.nansum(westward_transport_current * 80000))
        
    return Cell_net_transport_bw

def cell_mass_bilan(name ='Cell_A',area_lat = [75,77.5], area_lon = [-10,0]):
    """
        Compute the mass bilan of sea ice melted in the cell during each month by comparing the variation of sea ice along the month with the net export (or import) along the month
    """
    
    lon, lat, sit, sic, sit_uncertainty, time, area_marker = extracting_data_sit(area_lat = area_lat, area_lon = area_lon)
    
    SIV_variation_daily = []
    Cell_net_transport_bw = Net_transport_cell(name =name,area_lat = area_lat, area_lon = area_lon) #[m^3] positive when net import of sea ice on the cell
    Sea_ice_melt_on_the_cell = []
    SIV_daily = []

    useful_index = []
    day_of_year_with_sit_data = []
    day_gap = [] #Because we don't have daily data, we record in this dataset the number of day between the last element and the current element
    i = 0
    starting_date = 734419
    for time_ in time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
        if corresponding_date.year >= year_ and corresponding_date.year < year_end:
            if len(useful_index) == 0:
                    #Before taking the first data day of the year we take the last data day of the previous year in order to take the net SIV change for the first data day of the year.
                    useful_index.append(i-1)
            useful_index.append(i)
            
            day_of_year_with_sit_data.append(corresponding_date.day)
            day_gap.append(int(time[i])-int(time[i-1]))
        i += 1
    # All the sit data_array for the period of interest are merged in the following array
    recorded_sit = np.array([np.where(area_marker == 1, sit[n] * sic[n],np.nan) for n in useful_index])
    for day in range(1,len(recorded_sit)):
        SIV_variation = (np.nansum(recorded_sit[day])-np.nansum(recorded_sit[day-1])) * 80000**2 #[m^3] positive when sea ice volume increase over the cell
        SIV_variation_per_day = SIV_variation/day_gap[day-1]

        Sea_ice_melt_on_the_cell.append((Cell_net_transport_bw[day-1] - SIV_variation_per_day)*1e-9) #[km^3] positive when sea ice melt, negative when sea ice form (could be seen as fresh water flux through water)
        #print(Sea_ice_melt_on_the_cell[-1][-1])
        SIV_variation_daily.append(SIV_variation_per_day * 1e-9)#[km^3]
        SIV_daily.append(np.nansum(recorded_sit[day]*80000**2 * 1e-9)) #[km^3]
    
    """ plt.subplot(311)
    plt.plot(range(len(SIV_daily)),SIV_daily)
    plt.title('siv')
    plt.grid()
    plt.subplot(312)
    plt.plot(range(len(Sea_ice_melt_on_the_cell)),Sea_ice_melt_on_the_cell)
    plt.grid()
    plt.title('melting')
    plt.subplot(313)
    plt.grid()
    plt.plot(range(len(Cell_net_transport_bw)),Cell_net_transport_bw)
    plt.title('transport')
    plt.show()
     """

            
    np.savetxt('Data/bw/'+name+'/Cell_'+name[-1]+'_fw_daily.txt',Sea_ice_melt_on_the_cell)
    np.savetxt('Data/bw/'+name+'/Cell_'+name[-1]+'_siv_daily.txt',SIV_daily)
    np.savetxt('Data/bw/'+name+'/Cell_'+name[-1]+'_transp_daily.txt',Cell_net_transport_bw)
    return Sea_ice_melt_on_the_cell
year_, year_end = 2011,2021
cell_mass_bilan(name = 'Cell_A',area_lat=[77.5,80],area_lon = [-10,0])
cell_mass_bilan(name = 'Cell_B',area_lat=[75,77.5],area_lon = [-10,0])
cell_mass_bilan(name = 'Cell_C',area_lat=[72.5,75],area_lon = [-20,-10])
cell_mass_bilan(name = 'Cell_D',area_lat=[70,72.5],area_lon = [-20,-10])

""" Sea_ice_melt_on_the_cell = cell_mass_bilan()
print(np.shape(Sea_ice_melt_on_the_cell))
plt.plot(np.arange(len(Sea_ice_melt_on_the_cell[0])),Sea_ice_melt_on_the_cell[0])
plt.grid()
plt.show() """
