import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os 
from datetime import timedelta, date
from scipy import interpolate, stats
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import math

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
    v_ds.close
    u_ds.close

    print('##### - Data extracted - #####\n')
    return lon, lat, u_gos, v_gos, time

def interpole_and_save_ke():
    """
        Interpolate kinetic energy over the standard grid of cryosat and save it in Data/bw/KE
    """
    lon_gos, lat_gos, u_gos, v_gos,time_gos = extracting_data_gos()

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

    
    for day in range(len(recorded_ugos)):
        print(f'{day}/{len(recorded_ugos)}')
        points = []
        value_ugos = []
        value_vgos = []
        for line in range(len(lat_gos)):
            for col in range(len(lon_gos)):
                points.append([float(lat_gos[line]),float(lon_gos[col])])
                value_ugos.append(recorded_ugos[day][line,col])
                value_vgos.append(recorded_vgos[day][line,col])
        interp_ugos = interpolate.griddata(points,value_ugos,(lat,lon),method = 'linear')
        interp_vgos = interpolate.griddata(points,value_vgos,(lat,lon),method = 'linear')

        kinetic_energy = 1/2 * (interp_ugos**2 + interp_vgos**2)  #In [J/kg]
        print(np.shape(kinetic_energy))

        np.savetxt(f'Data/bw/KE/{int(date_data[day][0])}-{int(date_data[day][1])}-{int(date_data[day][2])}.txt',kinetic_energy)
        np.savetxt(f'Data/bw/u_gos/{int(date_data[day][0])}-{int(date_data[day][1])}-{int(date_data[day][2])}.txt',interp_ugos)
        np.savetxt(f'Data/bw/v_gos/{int(date_data[day][0])}-{int(date_data[day][1])}-{int(date_data[day][2])}.txt',interp_vgos)

def compute_correlation():
    """
        This function compute the correlation coefficient for each cells based on data stored in Data/bw
    """
    ke = []
    fw_flux = []
    siv_var = []
    transp = []
    for file in os.listdir('Data/bw/KE'):
        ke.append(np.loadtxt(f'Data/bw/KE/{file}')[1:-1,1:-1])
        fw_flux.append(np.loadtxt(f'Data/bw/FW_flux/{file}'))
        siv_var.append(np.loadtxt(f'Data/bw/SIV_var/{file}'))
        transp.append(np.loadtxt(f'Data/bw/Transport/{file}'))
    ke = np.nan_to_num(np.array(ke))
    fw_flux = np.nan_to_num(np.array(fw_flux))
    siv_var = np.nan_to_num(np.array(siv_var))
    transp = np.nan_to_num(np.array(transp))

    corr_ke_fw_flux = np.zeros(np.shape(ke[0]))
    p_value_ke_fw_flux = np.zeros(np.shape(ke[0]))

    corr_ke_siv_var = np.zeros(np.shape(ke[0]))
    p_value_ke_siv_var = np.zeros(np.shape(ke[0]))

    corr_ke_transp = np.zeros(np.shape(ke[0]))
    p_value_ke_transp = np.zeros(np.shape(ke[0]))
    for line in range(np.shape(corr_ke_fw_flux)[0]):
        for col in range(np.shape(corr_ke_fw_flux)[1]):
            corr_ke_fw_flux[line,col] = stats.pearsonr(ke[:,line,col],fw_flux[:,line,col])[0]
            p_value_ke_fw_flux[line,col] = stats.pearsonr(ke[:,line,col],fw_flux[:,line,col])[1]

            corr_ke_siv_var[line,col] = stats.pearsonr(ke[:,line,col],siv_var[:,line,col])[0]
            p_value_ke_siv_var[line,col] = stats.pearsonr(ke[:,line,col],siv_var[:,line,col])[1]

            corr_ke_transp[line,col] = stats.pearsonr(ke[:,line,col],transp[:,line,col])[0]
            p_value_ke_transp[line,col] = stats.pearsonr(ke[:,line,col],transp[:,line,col])[1]

    np.savetxt(f'Data/bw/r_ke_fw.txt',corr_ke_fw_flux)
    np.savetxt(f'Data/bw/p_ke_fw.txt',p_value_ke_fw_flux)

    np.savetxt(f'Data/bw/r_ke_siv_var.txt',corr_ke_siv_var)
    np.savetxt(f'Data/bw/p_ke_siv_var.txt',p_value_ke_siv_var)

    np.savetxt(f'Data/bw/r_ke_transp.txt',corr_ke_transp)
    np.savetxt(f'Data/bw/p_ke_transp.txt',p_value_ke_transp)

def compute_correlation_lag_time():
    """
        This function compute the correlation coefficient with different lag times spaced by 15 days (the time resolution) for each cells based on data stored in Data/bw
    """
    date_data = np.loadtxt('Data/bw\date.txt')
    for lag in range(1,5):
        ke = []
        fw_flux = []
        siv_var = []
        transp = []
        for file in os.listdir('Data/bw/KE'):
            ke.append(np.loadtxt(f'Data/bw/KE/{file}')[1:-1,1:-1])
            fw_flux.append(np.loadtxt(f'Data/bw/FW_flux/{file}'))
            siv_var.append(np.loadtxt(f'Data/bw/SIV_var/{file}'))
            transp.append(np.loadtxt(f'Data/bw/Transport/{file}'))
        ke = np.nan_to_num(np.array(ke))
        fw_flux = np.nan_to_num(np.array(fw_flux))
        siv_var = np.nan_to_num(np.array(siv_var))
        transp = np.nan_to_num(np.array(transp))

        #Lag time

        ke = ke[lag:,:,:]
        fw_flux = fw_flux[:-lag]
        siv_var = siv_var[:-lag]
        transp = transp[:-lag]
        corr_ke_fw_flux = np.zeros(np.shape(ke[0]))
        p_value_ke_fw_flux = np.zeros(np.shape(ke[0]))

        corr_ke_siv_var = np.zeros(np.shape(ke[0]))
        p_value_ke_siv_var = np.zeros(np.shape(ke[0]))

        corr_ke_transp = np.zeros(np.shape(ke[0]))
        p_value_ke_transp = np.zeros(np.shape(ke[0]))

        for line in range(np.shape(corr_ke_fw_flux)[0]):
            for col in range(np.shape(corr_ke_fw_flux)[1]):
                corr_ke_fw_flux[line,col] = stats.pearsonr(ke[:,line,col],fw_flux[:,line,col])[0]
                p_value_ke_fw_flux[line,col] = stats.pearsonr(ke[:,line,col],fw_flux[:,line,col])[1]

                corr_ke_siv_var[line,col] = stats.pearsonr(ke[:,line,col],siv_var[:,line,col])[0]
                p_value_ke_siv_var[line,col] = stats.pearsonr(ke[:,line,col],siv_var[:,line,col])[1]

                corr_ke_transp[line,col] = stats.pearsonr(ke[:,line,col],transp[:,line,col])[0]
                p_value_ke_transp[line,col] = stats.pearsonr(ke[:,line,col],transp[:,line,col])[1]

        np.savetxt(f'Data/bw/lag/r_ke_fw_lag={date_data[lag,-1]}_days.txt',corr_ke_fw_flux)
        np.savetxt(f'Data/bw/lag/p_ke_fw_lag={date_data[lag,-1]}_days.txt',p_value_ke_fw_flux)

        np.savetxt(f'Data/bw/lag/r_ke_siv_var_lag={date_data[lag,-1]}_days.txt',corr_ke_siv_var)
        np.savetxt(f'Data/bw/lag/p_ke_siv_var_lag={date_data[lag,-1]}_days.txt',p_value_ke_siv_var)

        np.savetxt(f'Data/bw/lag/r_ke_transp_lag={date_data[lag,-1]}_days.txt',corr_ke_transp)
        np.savetxt(f'Data/bw/lag/p_ke_transp_lag={date_data[lag,-1]}_days.txt',p_value_ke_transp)

def plot_r_p():
    """
        This function plots and save correlation map based on data in Data/bw/r_... and Data/bw/p_...
    """
    lon, lat, sit, sic, time = extracting_data_sit()
    ke = []
    for file in os.listdir('Data/bw/KE'):
        ke.append(np.loadtxt('Data/bw/KE/' + file)[1:-1,1:-1])
    ke = np.array(ke)
    mean_ke = np.nanmean(ke,axis = 0)
    ######## - KE vs FW_flux - #########
    r = np.loadtxt('Data/bw/r_ke_fw.txt')
    p = np.loadtxt('Data/bw/p_ke_fw.txt')
    lon = lon[1:-1,1:-1]
    lat = lat[1:-1,1:-1]
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
    #fig, axs = plt.plots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
    #axs.set_extent([-47, 16, 60, 85], crs = ccrs.PlateCarree())
    
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    #cs_ = axs.contour(lon, lat, p,[0.01], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    cs__ = axs.contour(lon, lat, mean_ke,[0.5*np.nanmean(mean_ke)], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between kinetic energy and freshwater flux", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots/correlation/map/correlation_ke_fw.png")

    ######## - KE vs siv_variation - #########
    r = np.loadtxt('Data/bw/r_ke_siv_var.txt')
    p = np.loadtxt('Data/bw/p_ke_siv_var.txt')
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between kinetic energy and siv variation", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots/correlation/map/correlation_ke_siv_var.png")

    ######## - KE vs transport - #########
    r = np.loadtxt('Data/bw/r_ke_transp.txt')
    p = np.loadtxt('Data/bw/p_ke_transp.txt')
    fig = plt.figure(figsize=(9,7))
    axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -18))
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
    levels = np.linspace(-0.5,0.5,1000)
    cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
    axs.set_title("Correlation between kinetic energy and net sea ice transport", fontsize = 15)
    cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
    
    cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
    cb.ax.tick_params(labelsize=25)
    plt.savefig(f"Plots/correlation/map/correlation_ke_transp.png")

def plot_r_p_lag():
    """
        Plot and save the correlation maps for different lag_time
    """
    lon, lat, sit, sic, time = extracting_data_sit()
    date_data = np.loadtxt('Data/bw\date.txt')

    lon = lon[1:-1,1:-1]
    lat = lat[1:-1,1:-1]
    for lag in range(1,5):
        r = np.loadtxt(f'Data/bw/lag/r_ke_fw_lag={date_data[lag,-1]}_days.txt')
        p = np.loadtxt(f'Data/bw/lag/p_ke_fw_lag={date_data[lag,-1]}_days.txt')
        fig = plt.figure(figsize=(9,7))
        axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -15))
        #fig, axs = plt.plots(nrows = 1, ncols = 1, figsize=figsize, subplot_kw={'projection': projection})
        #axs.set_extent([-47, 16, 60, 85], crs = ccrs.PlateCarree())
        
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
        levels = np.linspace(-0.5,0.5,1000)
        cs = axs.contourf(lon, lat, r, levels = levels, cmap = "cmo.balance", transform=ccrs.PlateCarree())
        cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
        #cs_ = axs.contour(lon, lat, p,[0.01], colors = 'green',linestyles = 'dashed', transform=ccrs.PlateCarree())
        axs.set_title(f"Correlation between kinetic energy and freshwater flux\n with lag time = {int(date_data[lag,-1])} days", fontsize = 15)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = [-0.5,0,0.5])
        cb.ax.tick_params(labelsize=25)
        plt.savefig(f"Plots/correlation/map/lag/correlation_ke_fw_lag={date_data[lag,-1]}_days.png")

def plot_fw_flux():
    """
        For each month, plots the mean freshwater flux
    """
    lon, lat, sit, sic, time = extracting_data_sit()
    date_data = np.loadtxt('Data/bw/date.txt')
    previous = 10
    fw_flux_month = []
    fw_monthly_mean = []
    for i in range(len(date_data)):
        if previous == date_data[i,1]:
            fw_flux_month.append(np.loadtxt(f'Data/bw/FW_flux/{int(date_data[i,0])}-{int(date_data[i,1])}-{int(date_data[i,2])}.txt'))
            previous = date_data[i,1]

        else:
            fw_monthly_mean.append(np.nanmean(fw_flux_month,axis = 0))
            fw_flux_month = [np.loadtxt(f'Data/bw/FW_flux/{int(date_data[i,0])}-{int(date_data[i,1])}-{int(date_data[i,2])}.txt')]
            previous = date_data[i,1]

    lon = lon[1:-1,1:-1]
    lat = lat[1:-1,1:-1]
    fw_monthly_mean = np.array(fw_monthly_mean)
    fw_monthly_mean[fw_monthly_mean == 0] = np.nan
    for month in range(len(fw_monthly_mean)):
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
        fw_monthly_mean[month] = fw_monthly_mean[month]*1e-9
        min = -4
        max = 4
        levels = np.linspace(min,max,20)
        cs = axs.contourf(lon, lat, fw_monthly_mean[month],levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
        #cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
        
        current_month = int(np.mod(month+9,12)) +1
        current_year = int(2010 + math.floor((month+9)/12))
        axs.set_title(f"Freshwater flux {current_month}/{current_year}", fontsize = 15)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = np.arange(min,max+1))
        cb.ax.tick_params(labelsize=25)
        
        print(f'{current_month} - {current_year}')
        try:
            os.makedirs(f'Plots/mean/FW_flux/{current_year}')
        except:
            pass
        plt.savefig(f'Plots/mean/FW_flux/{current_year}/{current_year}-{current_month}.png')
        plt.close()

def plot_siv_variation():
    """
        For each month, plots the mean freshwater flux
    """
    lon, lat, sit, sic, time = extracting_data_sit()
    date_data = np.loadtxt('Data/bw/date.txt')
    previous = 10
    siv_var_month = []
    siv_var_monthly_mean = []
    for i in range(len(date_data)):
        if previous == date_data[i,1]:
            siv_var_month.append(np.loadtxt(f'Data/bw/SIV_var/{int(date_data[i,0])}-{int(date_data[i,1])}-{int(date_data[i,2])}.txt'))
            previous = date_data[i,1]

        else:
            siv_var_monthly_mean.append(np.nanmean(siv_var_month,axis = 0))
            siv_var_month = [np.loadtxt(f'Data/bw/FW_flux/{int(date_data[i,0])}-{int(date_data[i,1])}-{int(date_data[i,2])}.txt')]
            previous = date_data[i,1]

    lon = lon[1:-1,1:-1]
    lat = lat[1:-1,1:-1]
    siv_var_monthly_mean = np.array(siv_var_monthly_mean)
    siv_var_monthly_mean[siv_var_monthly_mean == 0] = np.nan
    for month in range(len(siv_var_monthly_mean)):
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
        siv_var_monthly_mean[month] = siv_var_monthly_mean[month]*1e-9
        min = -4
        max = 4
        levels = np.linspace(min,max,20)
        cs = axs.contourf(lon, lat, siv_var_monthly_mean[month],levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
        #cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
        
        current_month = int(np.mod(month+9,12)) +1
        current_year = int(2010 + math.floor((month+9)/12))
        axs.set_title(f"SIV variation {current_month}/{current_year}", fontsize = 15)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = np.arange(min,max+1))
        cb.ax.tick_params(labelsize=25)
        
        print(f'{current_month} - {current_year}')
        try:
            os.makedirs(f'Plots/mean/SIV_var/{current_year}')
        except:
            pass
        plt.savefig(f'Plots/mean/SIV_var/{current_year}/{current_year}-{current_month}.png')
        plt.close()

def plot_transp():
    """
        For each month, plots the mean freshwater flux
    """
    lon, lat, sit, sic, time = extracting_data_sit()
    date_data = np.loadtxt('Data/bw/date.txt')
    previous = 10
    transp_month = []
    transp_monthly_mean = []
    for i in range(len(date_data)):
        if previous == date_data[i,1]:
            transp_month.append(np.loadtxt(f'Data/bw/Transport/{int(date_data[i,0])}-{int(date_data[i,1])}-{int(date_data[i,2])}.txt'))
            previous = date_data[i,1]

        else:
            transp_monthly_mean.append(np.nanmean(transp_month,axis = 0))
            transp_month = [np.loadtxt(f'Data/bw/Transport/{int(date_data[i,0])}-{int(date_data[i,1])}-{int(date_data[i,2])}.txt')]
            previous = date_data[i,1]

    lon = lon[1:-1,1:-1]
    lat = lat[1:-1,1:-1]
    transp_monthly_mean = np.array(transp_monthly_mean)
    transp_monthly_mean[transp_monthly_mean == 0] = np.nan
    for month in range(len(transp_monthly_mean)):
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
        transp_monthly_mean[month] = transp_monthly_mean[month]*1e-9
        min = -4
        max = 4
        levels = np.linspace(min,max,20)
        cs = axs.contourf(lon, lat, transp_monthly_mean[month],levels = levels,extend = 'both',cmap = "cmo.balance", transform=ccrs.PlateCarree())
        #cs_ = axs.contour(lon, lat, p,[0.05], colors = 'red', transform=ccrs.PlateCarree())
        
        current_month = int(np.mod(month+9,12)) +1
        current_year = int(2010 + math.floor((month+9)/12))
        axs.set_title(f"Sea ice transport {current_month}/{current_year}", fontsize = 15)
        cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
        
        cb = plt.colorbar(cs, cax = cax, ticks = np.arange(min,max+1))
        cb.ax.tick_params(labelsize=25)
        
        print(f'{current_month} - {current_year}')
        try:
            os.makedirs(f'Plots/mean/Transport/{current_year}')
        except:
            pass
        plt.savefig(f'Plots/mean/Transport/{current_year}/{current_year}-{current_month}.png')
        plt.close()

if True:   
    year_,year_end = 2010,2021
    lon, lat, sit, sic, time = extracting_data_sit()
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
plot_transp()
plot_siv_variation()