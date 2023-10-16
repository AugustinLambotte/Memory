import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
"""
    This script is used to plot on a map the value of the trend over all the time span for the EGC energy and the fresh water flux and for each pixel
"""

def trend():
    """
        return maps of the trend of fw_flux and EGC energy
    """
    lon = np.loadtxt('Data/bw/lon.txt')[1:-1,1:-1]
    lat = np.loadtxt('Data/bw/lat.txt')[1:-1,1:-1]
    
    fw_flux = []
    ke = []
    sit = []
    sic = []
    X_drift = []
    Y_drift = []
    date_data = np.loadtxt('Data/bw/date.txt')
    date_data_2011_to_2019 = []
    for date in date_data:
        if date[0] >= 2011 and date[0] < 2020:
            date_data_2011_to_2019.append(date)
    date_data_2011_to_2019 = np.array(date_data_2011_to_2019)
    number_of_days = date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1]
    for date in date_data_2011_to_2019:
        fw_flux.append(np.loadtxt(f'Data/bw/FW_flux/{int(date[0])}-{int(date[1])}-{int(date[2])}.txt'))
        ke.append(np.loadtxt(f'Data/bw/KE/{int(date[0])}-{int(date[1])}-{int(date[2])}.txt')[1:-1,1:-1])
        sit.append(np.loadtxt(f'Data/bw/sit/{int(date[0])}-{int(date[1])}-{int(date[2])}.txt')[1:-1,1:-1])
        sic.append(np.loadtxt(f'Data/bw/sic/{int(date[0])}-{int(date[1])}-{int(date[2])}.txt')[1:-1,1:-1])
        X_drift.append(np.loadtxt(f'Data/bw/X_drift/{int(date[0])}-{int(date[1]):02d}-{int(date[2]):02d}.txt')[1:-1,1:-1])
        Y_drift.append(np.loadtxt(f'Data/bw/Y_drift/{int(date[0])}-{int(date[1]):02d}-{int(date[2]):02d}.txt')[1:-1,1:-1])

    fw_flux = np.array(fw_flux)
    ke = np.nan_to_num(np.array(ke))
    mean_ke = np.nanmean(ke,axis = 0)
    sic= np.nan_to_num(np.array(sic))
    sit= np.nan_to_num(np.array(sit))
    X_drift= np.nan_to_num(np.array(X_drift))
    Y_drift= np.nan_to_num(np.array(Y_drift))
    
    trend_fw = np.zeros((np.shape(fw_flux)[1:]))
    trend_ke = np.zeros((np.shape(fw_flux)[1:]))
    trend_sit = np.zeros((np.shape(fw_flux)[1:]))
    trend_sic = np.zeros((np.shape(fw_flux)[1:]))
    trend_X_drift = np.zeros((np.shape(fw_flux)[1:]))
    trend_Y_drift = np.zeros((np.shape(fw_flux)[1:]))
    for l in range(np.shape(trend_fw)[0]):
        for c in range(np.shape(trend_fw)[1]):
            slope = np.polyfit(date_data_2011_to_2019[:,-1],fw_flux[:,l,c], deg = 1)[0]
            trend_fw[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = np.polyfit(date_data_2011_to_2019[:,-1],ke[:,l,c], deg = 1)[0]
            trend_ke[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])
            
            slope = np.polyfit(date_data_2011_to_2019[:,-1],sit[:,l,c], deg = 1)[0]
            trend_sit[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])
            
            slope = np.polyfit(date_data_2011_to_2019[:,-1],sic[:,l,c], deg = 1)[0]
            trend_sic[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = np.polyfit(date_data_2011_to_2019[:,-1],X_drift[:,l,c], deg = 1)[0]
            trend_X_drift[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

            slope = np.polyfit(date_data_2011_to_2019[:,-1],Y_drift[:,l,c], deg = 1)[0]
            trend_Y_drift[l,c] = slope * (date_data_2011_to_2019[-1,-1] - date_data_2011_to_2019[0,-1])

    trend_fw *= 1e-9
    trend_ke[trend_ke == 0] = np.nan
    trend_fw[trend_fw == 0] = np.nan
    trend_sit[trend_sit == 0] = np.nan
    trend_sic[trend_sic == 0] = np.nan
    trend_X_drift[trend_X_drift == 0] = np.nan
    trend_Y_drift[trend_Y_drift == 0] = np.nan
    mean_ke[mean_ke == 0] = np.nan

    trend_drift = trend_X_drift**2 + trend_Y_drift**2
    return trend_fw, trend_ke, trend_sit, trend_sic, trend_drift, mean_ke, lon, lat 

def plot():
    """
        plot the results of trend
    """
    trend_fw, trend_ke, trend_sit, trend_sic, trend_drift, mean_ke, lon, lat = trend()

    #### - Freshwater flux trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
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
    levels = np.linspace(-1,1,100)
    cs = axs.contourf(lon, lat, trend_fw, levels=levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    
    axs.set_title(f"Trend of fresh water flux from 2011 to 2019", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1,-0.5,0,0.5,1])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/trend/FW_flux.png')
    plt.close()

    #### - kinetic energy trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
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

    trend_ke *= 1e+3 #Passing from J/kg to 1e-3J/kg
    levels = np.linspace(-10,10,100)
    cs = axs.contourf(lon, lat, trend_ke, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    
    axs.set_title(f"Trend of kinetic energy [mJ/kg] from 2011 to 2019", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-10,-5,0,5,10])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/trend/KE.png')
    plt.close()

    #### - SIT trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
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

    levels = np.linspace(-1.2,1.2,100)
    cs = axs.contourf(lon, lat, trend_sit, levels=levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    
    axs.set_title(f"Trend of SIT [m] from 2011 to 2019", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-1.2,-1,-0.5,0,0.5,1,1.2])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/trend/SIT.png')
    plt.close()
    
    #### - SIC trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
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
    levels = np.linspace(-0.35,0.35,100)
    cs = axs.contourf(lon, lat, trend_sic, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    
    axs.set_title(f"Trend of SIC from 2011 to 2019", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.3,-0.2,-0.1,0,0.1,0.2,0.3])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/trend/SIC.png')
    plt.close()

    #### - SID trend - ####

    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    
    xlim = [-35, 12]
    ylim = [65, 81]
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
    levels = np.linspace(-40,40,200)
    cs = axs.contourf(lon, lat, trend_drift, levels = levels, extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
    
    axs.set_title(f"Trend of SID from 2011 to 2019", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-40,-30,0,30,40])
    cb.ax.tick_params(labelsize=20)
    
    plt.savefig(f'Plots/trend/SID.png')
    plt.close()
plot()