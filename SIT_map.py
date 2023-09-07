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
import matplotlib.path as mpath
"""import rioxarray as riox
from shapely.geometry import Polygon"""

####### - Functions definitions - #####

def plot_sea_ice_map_september_march( projection, figsize=(9, 7), nrows = 4, ncols = 4):
    """ 
        Return plot of sea ice thickness for all february and march month.
    """

    date = [[2011,9,15],[2012,3,15],[2012,9,15],[2013,3,15],[2013,9,15],[2014,3,15],[2014,9,15],[2015,3,15],[2015,9,15],[2016,3,15]
            ,[2016,9,15],[2017,3,15],[2017,9,15],[2018,3,15],[2018,9,15],[2019,3,15]]
    fig, axs = plt.subplots(nrows = nrows, ncols = ncols, figsize=figsize, subplot_kw={'projection': projection})
    cs_mem = [] #To have only one colorbar for lal plot we have to put in the argument of fig.colorbar() the cs with maximum
                #ice thickness so we record the cs in cs_mem to call the cs of 2012-3-15 which has the maximum thickness value.
    for ax, date_plot in zip(axs.flat, date):
        ax.set_extent([-40, 10, 60, 85], crs = ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines()
        cs = ax.contourf(lon, lat, sit[give_time_snapshot(date_plot[0],date_plot[1],date_plot[2])], vmin = 0, vmax= 8, transform=ccrs.PlateCarree())
        cs_mem.append(cs)
        ax.set_title("{}-{}".format(date_plot[0],date_plot[1]))
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.suptitle('SIT [m]')
    fig.colorbar(cs_mem[1], ax=axs.ravel().tolist(), shrink=0.95)

    plt.show()

def plot_sea_ice_map_single_date(year, month, day, projection, figsize = (9,7)):
    """
        Plot a single plot of the sea_ice_thickness at the given date value.
    """
    time_snapshot = give_time_snapshot(year,month,day)
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
    axs.set_extent([-70, 25, 55, 85], crs = ccrs.PlateCarree())
    axs.coastlines()
    axs.gridlines()
    cs = axs.contourf(lon, lat, sit[time_snapshot], cmap='cmo.ice', transform=ccrs.PlateCarree())
    axs.set_title("{}/{}/{}".format(year,month,day))
    fig.colorbar(cs, ax = axs)
    plt.show()

def give_date(time_snapshot): 
    """  
        Return the date given a number n in time[n] in [0,235]
        time_coverage_start = "2010-10-01"
        time_coverage_end = "2020-07-31" 
    """
    time_coverage = time[-1] - time[0] #In days
    time_step = time_coverage/len(time)
    days_from_start = time_step * time_snapshot

    return date(2010,10,1) + timedelta(days=int(days_from_start))

def give_time_snapshot(year,month,day):
    """   
        Return the n indix of time[n] given a date. Bug to fix : give_date(give_time_snapshot(date)) doesn't return exactly the date.
    """
    days_from_start = (date(year,month,day) - date(2010,10,1)).days
    time_coverage = time[-1] - time[0] # In days
    time_step = time_coverage/len(time) # In days
    time_snapshot = days_from_start/time_step
    return int(time_snapshot)

def extracting_data(file = "C:/Users/Augustin/Downloads/ubristol_cryosat2_seaicethickness_nh_80km_v1p7.nc", lat_range = [64,80], lon_range = [-40,20]):
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
    lon = ds['Longitude']
    lat = ds['Latitude']
    sit = ds['Sea_Ice_Thickness'].where((ds.Sea_Ice_Thickness != 0) & (ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))
    sic = ds['Sea_Ice_Concentration'].where((ds.Sea_Ice_Thickness != 0) & (ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))
    sit_uncertainty = ds['Sea_Ice_Thickness_Uncertainty'].where((ds.Longitude > lon_min) & (ds.Longitude < lon_max) & (ds.Latitude > lat_min) & (ds.Latitude < lat_max))
    time =  ds['Time']
    print(ds.attrs)
    ds.close
    return lon, lat, sit, sic, sit_uncertainty, time

def sea_ice_volume_march_sept():
    """ 
        Sea_ice_volume plotted each year during september and march.
    """
    #date_ = [date(y,month,15) for y in range(2011,2019)]
    march_SIV_per_year = []
    september_SIV_per_year = []
    surface_grid_cell = 80000**2 #resolution of 80km.
    for year in range(2011,2019):
        march_SIV_per_year.append(np.nansum(mensual_mean(year,3)) * surface_grid_cell)
        september_SIV_per_year.append(np.nansum(mensual_mean(year,9)) * surface_grid_cell)

    plt.plot([x for x in range(2011,2019)], march_SIV_per_year, label = "March")
    plt.plot([x for x in range(2011,2019)], september_SIV_per_year, label = "September")
    plt.legend()
    plt.xlabel('year')
    plt.ylabel('Sea_ice_volume [m^3]')
    plt.grid()
    plt.title("Mean SIV in March and September [m^3]")
    plt.show()

def mensual_mean(year, month):
    """
        Return the mensual mean value of SIT
    """
    useful_date = []
    useful_index = []
    starting_date = 734419
    i = 0
    for time_ in time:
        corresponding_date = date(2010,10,1) + timedelta(days = int(time_) - starting_date)
        if corresponding_date.year == year and corresponding_date.month == month:
            useful_date.append(int(time_))
            useful_index.append(i)
            print(corresponding_date) 
        i += 1
    recorded_sit = np.array([sit.sel(t = n) * sic.sel(t = n) for n in useful_index])
    mean_sit = np.nanmean(recorded_sit, axis = 0) 
    return mean_sit

def plot_mensual_mean(year, month, projection, figsize = (9,7), save = False):
    """
        Plot a single plot of the sea_ice_thickness at the given date value.
    """
    if save:
        dates = [[2011,3,15],[2011,9,15],[2012,3,15],[2012,9,15],[2013,3,15],[2013,9,15],[2014,3,15],[2014,9,15],[2015,3,15],[2015,9,15],[2016,3,15]
                ,[2016,9,15],[2017,3,15],[2017,9,15],[2018,3,15],[2018,9,15],[2019,3,15],[2019,9,15]]
        for year in range(2011,2020):
            for month in range(1,13):
                date = [year,month]
                print(f"### - Saving SIT: {date[0]}-{date[1]} - ###\n")
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
                cs = axs.contourf(lon, lat, mensual_mean(date[0],date[1]), levels = levels, cmap = "cmo.ice", transform=ccrs.PlateCarree())
                axs.set_title("Monthly averaged SIT value in meters for {}/{}".format(date[0],date[1]))
                
                cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0,0.02,axs.get_position().height])
                fig.colorbar(cs, cax = cax, ticks = [0,1,2,3,4])
                plt.grid
                plt.savefig(f"Plots/mean/Sea_ice/{year}/SIT_mean_{date[0]}-{date[1]}.png")
    else:
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
        axs.set_extent([-70, 25, 55, 85], crs = ccrs.PlateCarree())
        axs.coastlines()
        axs.gridlines()

        levels = np.linspace(0,7,1000)
        cs = axs.contourf(lon, lat, mensual_mean(year,month), levels= levels, cmap = "cmo.ice", transform=ccrs.PlateCarree())
        axs.set_title("Mensual mean SIT value for {}/{}".format(year,month))

        fig.colorbar(cs, ax = axs)
        plt.show()
    
##### - Main - #####
lon, lat, sit, sic, sit_uncertainty, time = extracting_data()
# The projection of the map
""" myProj = ccrs.LambertConformal(central_longitude=-20)
myProj._threshold = myProj._threshold/20.  # Set for higher precision of the projection
projection = myProj """
projection = ccrs.LambertConformal(central_longitude = -18)
plot_mensual_mean(2015,5,projection,save = True)
