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

    
    #Keeping data only on the area of interest
    epsilon = 0
    area_lat = [75-epsilon,77.5+epsilon]
    area_lon = [-10-epsilon,0+epsilon]
    area_marker = np.where((lat>=area_lat[0]) & (lat<=area_lat[1]) & (lon>=area_lon[0]) & (lon<=area_lon[1]),1,0)
    
    lat = lat.where((area_marker == 1))
    sit = sit.where((area_marker == 1))
    sic = sic.where((area_marker == 1))
    lon = lon.where((area_marker == 1))
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time, area_marker

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
year_, year_end = 2011,2020
lon, lat, sit, sic, sit_uncertainty, time, area_marker = extracting_data_sit()
for year in range(year_,year_end):
    print(year)
    for month in month_names:
        for file in os.listdir(f'Data/X_drift/{year}/{month}'):
            X_drift[f'y{year}'][month].append(np.where(area_marker == 1,np.loadtxt(f'Data/X_drift/{year}/{month}/' + file),np.nan))
            Y_drift[f'y{year}'][month].append(np.where(area_marker == 1,np.loadtxt(f'Data/Y_drift/{year}/{month}/' + file),np.nan))
###


northward_transport = np.zeros((year_end -year_,12))
southward_transport = np.zeros((year_end -year_,12))
eastward_transport = np.zeros((year_end -year_,12))
westward_transport = np.zeros((year_end -year_,12))
Cell_mass_bilan = np.zeros((year_end -year_,12))
for year in range(year_,year_end):
    for month in range(1,13):
        month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]
        useful_date = []
        useful_index = []
        day_of_month_with_sit_data = []
        i = 0
        starting_date = 734419
        for time_ in time:
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

        northward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
        southward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
        eastward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))
        westward_si_drift = np.zeros(np.shape(recorded_si_drift_X[0]))

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
        for day in range(len(recorded_sit)):
            for line in range(len(recorded_sit[day])):
                for col in range(len(recorded_sit[day][0])):
                    if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line-1,col]): #in this case we are in the northerest side of the cell
                        #Minus sign because we ant thhe transport to be positive when it points into the cell
                        northward_si_drift[line,col] = -(recorded_si_drift_X[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360))
                    if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line,col-1]): #in this case we are in the easterest side of the cell
                        southward_si_drift[line,col] = recorded_si_drift_X[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360)
                    if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line,col+1]): #in this case we are in the westerest side of the cell
                        eastward_si_drift[line,col] = -(recorded_si_drift_X[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360))
                    if not np.isnan(recorded_sit[day][line,col]) and np.isnan(recorded_sit[day][line+1,col]): #in this case we are in the southerest side of the cell
                        #Minus sign because we ant thhe transport to be positive when it points into the cell
                        westward_si_drift[line,col] = recorded_si_drift_X[day][line,col] * np.cos(lon[line,col] * 2*np.pi/360) + recorded_si_drift_Y[day][line,col] * np.sin(lon[line,col] * 2*np.pi/360)
            
            #Transport_ are array covering the same surface as recorded_si_drift and recorded_si_drift. Where there is sea ice drift data, the cell is filled with siv*si_drift.
            northward_transport_current = np.where(abs(northward_si_drift) >0 , northward_si_drift * daily_interp_sit[day], np.nan)
            southward_transport_current = np.where(abs(southward_si_drift) >0 , southward_si_drift * daily_interp_sit[day], np.nan)
            eastward_transport_current = np.where(abs(eastward_si_drift) >0 , eastward_si_drift * daily_interp_sit[day], np.nan)
            westward_transport_current = np.where(abs(westward_si_drift) >0 , westward_si_drift * daily_interp_sit[day], np.nan)

            northward_transport[year-year_,month-1] += np.nansum(northward_transport_current * 80000 * 1000)
            southward_transport[year-year_,month-1] += np.nansum(southward_transport_current * 80000 * 1000)
            eastward_transport[year-year_,month-1] += np.nansum(eastward_transport_current * 80000 * 1000)
            westward_transport[year-year_,month-1] += np.nansum(westward_transport_current * 80000 * 1000)

            #In m^3/month
            Cell_mass_bilan[year-year_,month-1] = northward_transport[year-year_,month-1] + southward_transport[year-year_,month-1] + eastward_transport[year-year_,month-1] + westward_transport[year-year_,month-1]


np.savetxt('Data/Cell_A/Mass_bilan_A.txt',Cell_mass_bilan)

""" figsize = (10,5)
projection = ccrs.LambertConformal(central_longitude = -18)
fig = plt.figure(figsize=figsize)
axs = plt.axes(projection = projection)
#fig, axs = plt.plots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
#axs.set_extent([-47, 16, 60, 85], crs = ccrs.PlateCarree())

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
levels = np.linspace(0,4,1000)
mean_sit,mean_sic = sit[10],sic[10]
cs = axs.contourf(lon, lat, mean_sit, levels = levels, cmap = "cmo.ice", transform=ccrs.PlateCarree())
cs_ = axs.contour(lon, lat, mean_sic,[0.15], colors = 'red', transform=ccrs.PlateCarree())
cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])

cb = plt.colorbar(cs, cax = cax, ticks = [0,1,2,3,4])
cb.ax.tick_params(labelsize=25)
plt.show() """