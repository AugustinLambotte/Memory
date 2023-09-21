import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from datetime import datetime, date, timedelta
import matplotlib.path as mpath
import cmocean as cmo
import cartopy.crs as ccrs


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
        for month in range(1,2):
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
                dX_interp = interpolate.griddata(points, values_dX, (lat_sit, lon_sit), method='linear')
                dY_interp = interpolate.griddata(points, values_dY, (lat_sit, lon_sit), method='linear')

                
                X_drift[f'y{year}'][month_name].append(dX_interp)
                Y_drift[f'y{year}'][month_name].append(dY_interp)
    
    
    return X_drift, Y_drift, lat, lon



if __name__ == '__main__':
    
    year_ = 2011
    year_end = 2012
    #X_drift, Y_drift, lat, lon = extracting_SI_drift()
    lon_sit,lat_sit, sit, sic, sit_uncertainty, sit_time= extracting_data_sit()
    X_drift, Y_drift, lat_drift, lon_drift = extracting_SI_drift()
    X_sit,Y_sit = np.meshgrid(sit.x,sit.y)
    X_dri,Y_dri = np.meshgrid(lat_drift.xc,lat_drift.yc)
    northward_transport = np.zeros((year_end -year_,12))
    for year in range(year_,year_end):
        northward_transport_current_month = []
        for month in range(1,2):

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
            recorded_si_drift = np.array([X_drift[f"y{year}"][month_names[month-1]][day-1] for day in day_of_month_with_sit_data])
            
            
            for day in range(len(recorded_si_drift)):
                coord_siv = [] #Record the coord [lat,lon] of the points recorded in siv_mean
                coord_dX = [] #Record the coord [lat,lon] of the points recorded in dX_mean
                siv_fram = [] #This array record all the northest point with siv for each longitude.
                drift_fram = [] #This array record all the northest point with SI drift for each longitude.
                matrix_coord_dX = [] # record the coord in data coord [x,y] where there is the northest drift data. This is used to find 
                                    # which siv is over the cell  
                for y in sit.y:
                    for x in reversed(sit.x):
                        if not np.isnan(recorded_sit[day][y,x]):
                            siv_fram.append(recorded_sit[day][y,x])
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
                    northward_transport_current_day += drift_fram[i] * 1000 * recorded_sit[day][matrix_coord_dX[i][0],matrix_coord_dX[i][1]] * 80000 #Transport en [m^3/jour]
                northward_transport_current_month.append(northward_transport_current_day)
            northward_transport[year-year_,month-1] = np.mean(northward_transport_current_month)*30
        
        
        plt.figure(figsize=(10,6))
        plt.plot([month for month in range(1,13)],[northward_transport[year-year_,month]* 1e-9 for month in range(12)])
        plt.grid()
        plt.title(f"Northward transport of sea ice through Fram Strait during year {year} in [km^3]")
        plt.xlabel('month')
        plt.ylabel('[km^3]')
        plt.savefig(f"Plots/Fram_strait/2.0/{year}_Mass_bilan.png")

            
            
    """fig, axs = plt.subplots(nrows = 1, ncols = 2, subplot_kw={'projection': ccrs.LambertConformal(central_longitude = -18)})
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
    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs[0]) - axs[0].transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs[0].set_boundary(rect_in_target)
    axs[0].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
        
    axs[0].coastlines()
    axs[0].gridlines()
    levels = np.linspace(0,4,1000)
    cs = axs[0].contourf(lon_sit, lat_sit, mensual_mean_siv(year_,1), levels = levels, cmap = "cmo.ice", transform=ccrs.PlateCarree())
    axs[0].set_title("Monthly averaged SIT value in meters for {}".format(year))
        
    cax = fig.add_axes([axs[0].get_position().x1+0.01,axs[0].get_position().y0,0.02,axs[0].get_position().height])
    fig.colorbar(cs, cax = cax, ticks = [0,1,2,3,4])

    proj_to_data   = ccrs.PlateCarree()._as_mpl_transform(axs[1]) - axs[1].transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs[1].set_boundary(rect_in_target)
    axs[1].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])
        
    axs[1].coastlines()
    axs[1].gridlines()
    levels = np.linspace(0,0.4,1000)
    cs = axs[1].contourf(lon, lat, mensual_mean_drift(year_,1)[0], cmap = "cmo.ice", transform=ccrs.PlateCarree())
    axs[1].set_title("Monthly averaged SIT value in meters for {}".format(year_))
        
    cax = fig.add_axes([axs[1].get_position().x1+0.01,axs[1].get_position().y0,0.02,axs[1].get_position().height])
    fig.colorbar(cs, cax = cax, ticks = [0,0.1,0.2,0.3,0.4])
    plt.show()

    print(drift_fram)
    print(coord_dX)
    print('######')
    print(siv_fram)
    print(coord_siv)

    
    plt.subplot(211)
    plt.imshow(siv_mean)
    plt.title('siv')
    plt.subplot(212)
    plt.imshow(dX_mean)
    plt.title('drift')
    plt.show() """
            
    
    