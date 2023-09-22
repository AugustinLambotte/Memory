import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime, date, timedelta
import matplotlib.path as mpath
import cmocean as cmo
import cartopy.crs as ccrs
import scipy

def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [78,80.5], lon_range = [-40,20]):
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
    sit = ds['Sea_Ice_Thickness'].where((ds.Sea_Ice_Thickness != 0)) 
    sit = sit.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sic = ds['Sea_Ice_Concentration'].where((ds.Sea_Ice_Thickness != 0)) 
    sic = sic.where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max) & (ds.Latitude > 65.4 + (76.5-65.4)/(9+17) * (ds.Longitude + 17)), drop = True)
    time =  ds['Time']
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time

def extracting_SI_drift(lat_range = [78,80.5], lon_range = [-40,20]):
    
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
                    X_drift[f'y{year}'][month_name].append(dX_interp)
                    Y_drift[f'y{year}'][month_name].append(dY_interp)
                else:
                    print('Manual interpolation')
                    #In this case, there are not enough data points (less than 2) to perform a classical interpolation
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

                    X_drift[f'y{year}'][month_name].append(manualy_interpolated_X)
                    Y_drift[f'y{year}'][month_name].append(manualy_interpolated_Y)
        print(f"{day_ignored} days man interp during year {year}")
    
    return X_drift, Y_drift, lat, lon



if __name__ == '__main__':
    
    year_ = 2010
    year_end = 2021
    #X_drift, Y_drift, lat, lon = extracting_SI_drift()
    lon_sit,lat_sit, sit, sic, sit_uncertainty, sit_time= extracting_data_sit()
    X_drift, Y_drift, lat_drift, lon_drift = extracting_SI_drift()
    X_sit,Y_sit = np.meshgrid(sit.x,sit.y)
    X_dri,Y_dri = np.meshgrid(lat_drift.xc,lat_drift.yc)
    northward_transport = np.zeros((year_end -year_,12))
    for year in range(year_,year_end):
        northward_transport_current_month = []
        for month in range(1,13):

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
            recorded_sit = np.array([sit.sel(t = n) * sic.sel(t = n) for n in useful_index])
            recorded_si_drift = np.array(X_drift[f"y{year}"][month_names[month-1]][:]) 
            
            def interp_sit_daily():
                
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
                points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
                values_sit = []
                for day in day_of_month_with_sit_data:
                    for i in range(len(lat_sit)):
                        for j in range(len(lon_sit[0])):
                            if recorded_sit[day_of_month_with_sit_data.index(day)][i,j] !=0  and not np.isnan(recorded_sit[day_of_month_with_sit_data.index(day)][i,j]):
                                points.append([day,lat_sit[i,j],lon_sit[i,j]])
                                values_sit.append(recorded_sit[day_of_month_with_sit_data.index(day)][i,j])
                points = np.array(points)
                values_sit = np.array(values_sit)

                lat_3d = [lat_sit for day in range(1,nb_days+1)]
                lon_3d = [lon_sit for day in range(1,nb_days+1)]
                day_3d = [np.ones((len(lat_sit),len(lat_sit[0]))) * day for day in range(1,nb_days+1)]
                daily_interp_sit = interpolate.griddata(points, values_sit, (day_3d,lat_3d, lon_3d))
                return daily_interp_sit
            
            def interp_sit_daily_2():
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
                if len(day_of_month_with_sit_data) < 2:
                    return [recorded_sit[0] for day in range(nb_days)]
                else:
                    interp = interpolate.interp1d(day_of_month_with_sit_data,recorded_sit,axis = 0,fill_value="extrapolate")
                    return interp([day for day in range(1,nb_days+1)])
            
            daily_interp_sit = interp_sit_daily_2()
            for day in range(len(recorded_si_drift)):
                coord_siv = []       # Record the coord [lat,lon] of the points recorded in siv_mean
                coord_dX = []        # Record the coord [lat,lon] of the points recorded in dX_mean
                siv_fram = []        # This array record all the northest point with siv for each longitude.
                drift_fram = []      # This array record all the northest point with SI drift for each longitude.
                matrix_coord_dX = [] # record the coord in data coord [x,y] where there is the northest drift data. This is used to find 
                                     # which siv is over the cell  
                for y in sit.y:
                    for x in reversed(sit.x):
                        #print(daily_interp_sit[day][y,x])
                        if not np.isnan(daily_interp_sit[day][y,x]):
                            siv_fram.append(daily_interp_sit[day][y,x])
                            coord_siv.append([float(lat_sit[y,x]),float(lon_sit[y,x])])
                            break

                for y in sit.y:
                    for x in reversed(sit.x):
                        if not np.isnan(recorded_si_drift[day][y,x]):
                            drift_fram.append(recorded_si_drift[day][y,x])
                            coord_dX.append([float(lat_sit[y,x]),float(lon_sit[y,x])])
                            matrix_coord_dX.append([int(y),int(x)])
                            break

                northward_transport_current_day = 0
                for i in range(len(coord_dX)):
                    northward_transport_current_day += drift_fram[i] * 1000 * daily_interp_sit[day][matrix_coord_dX[i][0],matrix_coord_dX[i][1]] * 80000 #Transport en [m^3/jour]
                
                northward_transport[year-year_,month-1] += northward_transport_current_day
        
        
        """ plt.figure(figsize=(10,6))
        plt.plot([month for month in range(1,13)],[northward_transport[year-year_,month]* 1e-9 for month in range(12)],marker = 'v')
        plt.grid()
        plt.title(f"Northward transport of sea ice through Fram Strait during year {year} in [km^3]")
        plt.xlabel('month')
        plt.ylabel('[km^3]')
        plt.savefig(f"Plots/Fram_strait/3.0/{year}_Mass_bilan.png") """

    plt.figure(figsize=(10,6))
    plt.plot([year for year in range(year_,year_end)],[np.nansum(northward_transport[year-year_,:])* 1e-9 for year in range(year_,year_end)],marker = 'v')
    plt.grid()
    plt.title(f"Annual sea ice northward transport through Fram Strait in [km^3]")
    plt.xlabel('year')
    plt.ylabel('[km^3]')
    plt.savefig(f"Plots/Fram_strait/3.0/Annual_mean_mass_bilan.png")

            
            

            
    
    