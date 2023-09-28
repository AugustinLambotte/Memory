import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime, date, timedelta
import matplotlib.path as mpath
import cmocean as cmo
import cartopy.crs as ccrs
import scipy

def extracting_data_sit(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [78,82.5], lon_range = [-20,20]):
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
    """ plt.subplot(311)
    plt.imshow() """

    #Keeping data only on the gate
    gate_lat = 80
    gate_marker = np.zeros(lat.shape)
    for i in range(len(lat)):
        gate_marker[i,int((abs(lat[i,:]-gate_lat).argmin()))] = 1
    
    
    lat = lat.where((gate_marker == 1))
    sit = sit.where((gate_marker == 1))
    sic = sic.where((gate_marker == 1))
    lon = lon.where((gate_marker == 1))
    

    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time, gate_marker

def extracting_SI_drift(lat_range = [78,82.5], lon_range = [-20,20]):
    
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

                    #Keeping data only on the gate
                    dX_interp = np.where(gate_marker == 1, dX_interp,np.nan)
                    dY_interp = np.where(gate_marker == 1, dY_interp,np.nan)
                    
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

                    #Keeping data only on the gate
                    manualy_interpolated_X = np.where(gate_marker == 1, manualy_interpolated_X,np.nan)
                    manualy_interpolated_Y = np.where(gate_marker == 1, manualy_interpolated_Y,np.nan)

                    X_drift[f'y{year}'][month_name].append(manualy_interpolated_X)
                    Y_drift[f'y{year}'][month_name].append(manualy_interpolated_Y)
        print(f"{day_ignored} days man interp during year {year}")
    
    return X_drift, Y_drift, lat, lon



if __name__ == '__main__':
    
    year_ = 2011
    year_end = 2020
    #X_drift, Y_drift, lat, lon = extracting_SI_drift()
    lon_sit,lat_sit, sit, sic, sit_uncertainty, sit_time, gate_marker= extracting_data_sit()
    X_drift, Y_drift, lat_drift, lon_drift = extracting_SI_drift()
    X_sit,Y_sit = np.meshgrid(sit.x,sit.y)
    X_dri,Y_dri = np.meshgrid(lat_drift.xc,lat_drift.yc)
    northward_transport = np.zeros((year_end -year_,12))
    for year in range(year_,year_end):
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
            recorded_si_drift_X = np.array(X_drift[f"y{year}"][month_names[month-1]][:]) 
            recorded_si_drift_Y = np.array(Y_drift[f"y{year}"][month_names[month-1]][:]) 
               
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
                
                if len(day_of_month_with_sit_data) == 1:
                    return [recorded_sit[0] for day in range(nb_days)]
                else:
                    interp = interpolate.interp1d(day_of_month_with_sit_data,recorded_sit,axis = 0, bounds_error = False, fill_value = (recorded_sit[0],recorded_sit[-1]))
                    return interp([day for day in range(1,nb_days+1)])
            
            daily_interp_sit = interp_sit_daily()
            for day in range(len(recorded_si_drift_X)):

                northward_si_drift = recorded_si_drift_X[day] * np.sin(lon_sit * 2*np.pi/360) + recorded_si_drift_Y[day] * np.cos(lon_sit * 2*np.pi/360)

                #Transport is an array covering the same surface as recorded_si_drift and recorded_si_drift. Where there is sea ice drift data, the cell is filled with siv*si_drift.
                transport = np.where(abs(northward_si_drift) >0 , northward_si_drift * daily_interp_sit[day], np.nan)
                
                northward_transport[year-year_,month-1] += np.nansum(transport * 80000 * 1000)
        #Ricker et al. (2018) datas from 2011 to 2016
        Ricker = np.array([[-267,-21,-540,-279,np.nan,np.nan,np.nan,np.nan,np.nan,-164,-214,-354],
                  [-129,-381,-379,-487,np.nan,np.nan,np.nan,np.nan,np.nan,-203,-182,-187],
                  [-103,-163,-299,-318,np.nan,np.nan,np.nan,np.nan,np.nan,-215,-400,-231],
                  [-78,-195,-345,-452,np.nan,np.nan,np.nan,np.nan,np.nan,-200,-165,-373],
                  [-160,-425,-429,-354,np.nan,np.nan,np.nan,np.nan,np.nan,-52,-261,-275],
                  [-177,-352,-348,-310,np.nan,np.nan,np.nan,np.nan,np.nan,-129,-151,-307]])
        
        M1 = np.array([[-238,-24,-478,-255,np.nan,np.nan,np.nan,np.nan,np.nan,-149,-163,-293],
              [-109,-299,-287,-428,np.nan,np.nan,np.nan,np.nan,np.nan,-207,-157,-125],
              [-80,-122,-254,-254,np.nan,np.nan,np.nan,np.nan,np.nan,-212,-372,-211],
              [-49,-105,-240,-401,np.nan,np.nan,np.nan,np.nan,np.nan,-203,-122,-307],
              [-129,-358,-328,-284,np.nan,np.nan,np.nan,np.nan,np.nan,-72,-215,-243],
              [-129,-272,-255,-264,np.nan,np.nan,np.nan,np.nan,np.nan,-98,-90,-243]])
        
        M2 = np.array([[-238,-34,-442,-230,-278,-185,-115,-64,-28,-151,-175,-290],
              [-137,-300,-267,-372,-334,-218,-187,-131,-100,-160,-149,-136],
              [-78,-109,-217,-219,-194,-140,-107,-98,-26,-228,-367,-191],
              [-61,-114,-282,-425,-232,-161,-112,-184,-194,-170,-162,-283],
              [-129,-355,-339,-308,-171,-240,-114,-11,-107,-78,-192,-244],
              [-150,-267,-287,-289,-196,-194,-113,-198,-75,-97,-72,-222]])
        
        plt.figure(figsize=(12,7))
        plt.plot([month for month in range(1,13)],[northward_transport[year-year_,month]* 1e-9 for month in range(12)],marker = 'v', color = 'blue', label = 'myself')
        if year < 2017:
            plt.plot([month for month in range(1,13)],[Ricker[year-year_,month] for month in range(12)],marker = 'o', color = 'red',label = 'Ricker et al. (2018)')
            plt.plot([month for month in range(1,13)],[M1[year-year_,month] for month in range(12)],marker = 'o',color = 'green', label = 'Min et al.(2019) - M1')
            plt.plot([month for month in range(1,13)],[M2[year-year_,month] for month in range(12)],marker = 'o', color = 'orange', label = 'Min et al.(2019) - M2')
        plt.grid()
        plt.legend(fontsize = "15")
        plt.ylim(-550,20)
        plt.title(f"Northward transport of sea ice through Fram Strait during year {year} in [km^3]",fontdict ={'fontsize':21})
        plt.xlabel('month',fontdict = {'fontsize':20})
        plt.ylabel('[km^3]',fontdict = {'fontsize':20})
        plt.yticks(fontsize = 20)
        plt.xticks(fontsize = 20)
        plt.savefig(f"Plots/Fram_strait/5.0/Comparison/{year}_Mass_bilan.png")
    

    plt.figure(figsize=(12,7))
    plt.plot([year for year in range(year_,year_end)],[np.nansum(northward_transport[year-year_,:])* 1e-9 for year in range(year_,year_end)],marker = 'v',color = 'blue',label = 'myself')
    plt.plot([year for year in range(year_,2017)],[np.nansum(M2[year-year_,:]) for year in range(year_,2017)],marker = 'o', color = 'orange', label = 'Min et al.(2019) - M2')
    plt.legend(fontsize = "15")
    plt.grid()
    plt.title(f"Annual sea ice northward transport through Fram Strait in [km^3]",fontdict ={'fontsize':25})
    plt.xlabel('year',fontdict = {'fontsize':20})
    plt.ylabel('[km^3]',fontdict = {'fontsize':20})
    plt.yticks(fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.savefig(f"Plots/Fram_strait/5.0/Comparison/Annual_mean_mass_bilan.png")

    plt.figure(figsize=(12,7))
    plt.plot([year for year in range(year_,year_end)],[np.nansum(northward_transport[year-year_,4:9])* 1e-9 for year in range(year_,year_end)],marker = 'v',color = 'blue',label = 'myself')
    plt.plot([year for year in range(year_,2017)],[np.nansum(M2[year-year_,4:9]) for year in range(year_,2017)],marker = 'o', color = 'orange', label = 'Min et al.(2019) - M2')
    plt.legend(fontsize = "15")
    plt.grid()
    plt.title(f"Summer sea ice northward transport through Fram Strait in [km^3]",fontdict ={'fontsize':25})
    plt.xlabel('year',fontdict = {'fontsize':20})
    plt.ylabel('[km^3]',fontdict = {'fontsize':20})
    plt.yticks(fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.savefig(f"Plots/Fram_strait/5.0/Comparison/Summer_mean_mass_bilan.png")

    plt.figure(figsize=(12,7))
    plt.plot([year for year in range(year_,year_end)],[(np.nansum(northward_transport[year-year_,:4]) +  np.nansum(northward_transport[year-year_,9:]))* 1e-9 for year in range(year_,year_end)],marker = 'v',color = 'blue',label = 'myself')
    plt.plot([year for year in range(year_,2017)],[(np.nansum(M2[year-year_,:4]) + np.nansum(M2[year-year_,9:])) for year in range(year_,2017)],marker = 'o', color = 'orange', label = 'Min et al.(2019) - M2')
    plt.plot([year for year in range(year_,2017)],[np.nansum(Ricker[year-year_,:]) for year in range(year_,2017)],marker = 'o', color = 'red',label = 'Ricker et al. (2018)')
    plt.plot([year for year in range(year_,2017)],[np.nansum(M1[year-year_,:]) for year in range(year_,2017)],marker = 'o',color = 'green', label = 'Min et al.(2019) - M1')
    plt.legend(fontsize = "15")
    plt.grid()
    plt.title(f"Winter sea ice northward transport through Fram Strait in [km^3]",fontdict ={'fontsize':25})
    plt.xlabel('year',fontdict = {'fontsize':20})
    plt.ylabel('[km^3]',fontdict = {'fontsize':20})
    plt.yticks(fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.savefig(f"Plots/Fram_strait/5.0/Comparison/winter_mean_mass_bilan.png")
            
            

            
    
    