import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os 
from datetime import timedelta, date
from scipy import interpolate
"""
    This script compute the correlation coefficient for every pixels (80km^2) over 
    the EGC between the intensity of the EGC (KE) and the fresh water flux.
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
    #sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    time =  ds['Time']

    
    ds.close
    return lon, lat, sit, sic, time

def save_siv_transp_fw():
    # For each cell NOT on the edge of the grid, we compute the net transport (positive when net import of sea ice) by considering the four cells 
    # on the edges of the cell. We stock this in Cell_transport. In Cell_siv_var we stock the variation of sea ice volume for each cell. We divide
    # them by the number of days between the two date to have a 'daily' siv variation.

    Cell_transport = np.zeros((np.shape(X_drift)[0],np.shape(X_drift)[1]-2,np.shape(X_drift)[2]-2))
    Cell_siv_var = np.zeros((np.shape(X_drift)[0],np.shape(X_drift)[1]-2,np.shape(X_drift)[2]-2))
    for day in range(np.shape(X_drift)[0]):
        print(f'{int(date_data[day][0])}-{int(date_data[day][1])}-{int(date_data[day][2])}')
        for line in range(1,np.shape(X_drift)[1]-1):
            for col in range(1,np.shape(X_drift)[2]-1):

                northward_si_drift = X_drift[day][line+1,col] * np.sin(lon[line+1,col] * 2*np.pi/360) - Y_drift[day][line+1,col] * np.cos(lon[line+1,col] * 2*np.pi/360)
                westward_si_drift = X_drift[day][line,col-1] * np.cos(lon[line,col-1] * 2*np.pi/360) + Y_drift[day][line,col-1] * np.sin(lon[line,col-1] * 2*np.pi/360)
                eastward_si_drift = -(X_drift[day][line,col+1] * np.cos(lon[line,col+1] * 2*np.pi/360) + Y_drift[day][line,col+1] * np.sin(lon[line,col+1] * 2*np.pi/360))
                southward_si_drift = - X_drift[day][line-1,col] * np.sin(lon[line-1,col] * 2*np.pi/360) + Y_drift[day][line-1,col] * np.cos(lon[line-1,col] * 2*np.pi/360)

                northward_si_drift = np.nan_to_num(northward_si_drift)
                westward_si_drift = np.nan_to_num(westward_si_drift)
                eastward_si_drift = np.nan_to_num(eastward_si_drift)
                southward_si_drift = np.nan_to_num(southward_si_drift)

                northward_transport = northward_si_drift * 1000 * recorded_sit[day+1][line+1,col] * 80000 #[m^3] we take [day+1] because recorded_sit is one longer than the other array because recorded_sit[0] is the day before the first day of the time range oof interest.
                southward_transport = southward_si_drift * 1000 * recorded_sit[day+1][line-1,col] * 80000 #[m^3]
                eastward_transport = eastward_si_drift * 1000 * recorded_sit[day+1][line,col+1] * 80000 #[m^3]
                westward_transport = westward_si_drift * 1000 * recorded_sit[day+1][line,col-1] * 80000 #[m^3]            

                Cell_transport[day][line-1,col-1] = northward_transport + southward_transport + eastward_transport + westward_transport
                Cell_siv_var[day][line-1,col-1] = (recorded_sit[day+1][line,col] - recorded_sit[day][line,col])*80000**2/time_gap[day] #Positive when siv increase over the cell
        
    fresh_water_flux = Cell_transport - Cell_siv_var
    for date_,i in zip(date_data,range(len(date_data))):
        np.savetxt(f'Data/bw/FW_flux/{int(date_[0])}-{int(date_[1])}-{int(date_[2])}.txt',fresh_water_flux[i])
        np.savetxt(f'Data/bw/SIV_var/{int(date_[0])}-{int(date_[1])}-{int(date_[2])}.txt',Cell_siv_var[i])
        np.savetxt(f'Data/bw/Transport/{int(date_[0])}-{int(date_[1])}-{int(date_[2])}.txt',Cell_transport[i])

def extracting_data_gos(area_lat = [60,83], area_lon = [-40,20]):
    """ 
        Given the file path "C:/.../..." to the .nc file it returns the longitude, latitude, sea_ice_thickness and time 
        restricted on the area defined by lat_range and lon_range
    """
    print("\n##############################")
    print("##### - Extracting data - ####")
    print("##############################\n")

    lon_min = area_lon[0]
    lon_max = area_lon[1]
    lat_min = area_lat[0]
    lat_max = area_lat[1]

    v_file = "C:/Users/Augustin/Downloads/Northward_sea_water_velocity"
    u_file = "C:/Users/Augustin/Downloads/Eastward_sea_water_velocity"
    
    v_ds = xr.open_dataset(v_file, decode_times = False)
    u_ds = xr.open_dataset(u_file, decode_times = False)
    
    u_gos= u_ds['uo'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude < lat_max),drop = True)
    v_gos= v_ds['vo'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude < lat_max),drop = True)
    
    u_gos = u_gos.sel(depth = 0.494025)
    v_gos = v_gos.sel(depth = 0.494025)
    
    lon = v_ds['longitude']
    lat = v_ds['latitude']
    time =  v_ds['time']
    print(u_ds.attrs)
    v_ds.close
    u_ds.close

    print('##### - Data extracted - #####\n')
    return lon, lat, u_gos, v_gos, time

def interpole_and_save_ke():
    """
        Interpolate kinetic energy over the standard grid of cryosat and save it in Data/bw/KE
    """
    lon_gos, lat_gos, u_gos,v_gos,time_gos = extracting_data_gos()

    #The following lines are selecting the date when there are sit data
    starting_date = 734419
    rec_date_sit = []
    for time_ in time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
        if corresponding_date.year >= year_ and corresponding_date.year < year_end:
            rec_date_sit.append(corresponding_date)

    #The following lines are selecting the corresponding days in the gos current data
    useful_index = []
    i = 0
    for time_ in time_gos:
        corresponding_date = date(1950,1,1) + timedelta(hours = int(time_)) 
        if corresponding_date in rec_date_sit:
            useful_index.append(i)
        i+=1

    # Creation of arrays with only the gos data when there is cryosat data
    recorded_ugos = np.array([u_gos.isel(time = n) for n in useful_index])
    recorded_vgos = np.array([v_gos.isel(time = n) for n in useful_index])

    #Interpolation over the cryosat spatial grid
    interp_ugos = []
    interp_vgos = []
    for day in range(len(recorded_ugos)):
        points = []
        value_ugos = []
        value_vgos = []
        for line in range(len(recorded_ugos[0])):
            for col in range(len(recorded_ugos[0][0])):
                points.append([lat_gos[line,col],lon_gos[line,col]])
                value_ugos.append(recorded_ugos[day][line,col])
                value_vgos.append(recorded_vgos[day][line,col])
        interp_ugos.append(interpolate.griddata(points,value_ugos,(lat,lon),method = 'linear'))
        interp_vgos.append(interpolate.griddata(points,value_vgos,(lat,lon),method = 'linear'))
    print(np.shape(interp_ugos))
    print(np.shape(interp_vgos))
    kinetic_energy = [np.nanmean(1/2 * (recorded_ugos[day]**2 + recorded_vgos[day]**2)) for day in range(len(recorded_ugos))] #In [J/kg]
    print(len(kinetic_energy))

    np.savetxt("Data/bw/"+name+"/Cell_"+name[-1]+"_index_daily.txt",kinetic_energy) 

year_,year_end = 2011,2021
lon, lat, sit, sic, time = extracting_data_sit()
interpole_and_save_ke()
X_drift = []
Y_drift = []
date_data = np.loadtxt('Data/bw\date.txt')
for file in os.listdir(f'Data/bw/X_drift'):
    X_drift.append(np.loadtxt(f'Data/bw/X_drift/{file}'))
    Y_drift.append(np.loadtxt(f'Data/bw/Y_drift/{file}'))


useful_index = []
day_of_year_with_sit_data = []
i = 0
starting_date = 734419
time_gap = []#record the time gap in days between each data points. (14, 15 or 16)

for time_ in time:
    corresponding_date = date(2010,10,1) + timedelta(days = int(time_)-starting_date)
    if corresponding_date.year >= year_ and corresponding_date.year < year_end:
        if len(useful_index) == 0:#We record one day before the starting date to be able to compute the SIV_variation for the first day.
            useful_index.append(i-1)
            
        useful_index.append(i)
        time_gap.append(int(time[i]) - int(time[i-1]))
        day_of_year_with_sit_data.append(corresponding_date.day)
    i += 1

# All the sit data_array for the month of interest are merged in the following array
recorded_sit = np.nan_to_num(np.array([sit[n] * sic[n] for n in useful_index]))
