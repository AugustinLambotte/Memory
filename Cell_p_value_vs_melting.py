import numpy as np
import matplotlib.pyplot as plt

Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']


#Total fresh_water_flux
for cell in Cell_name:
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('r',color = 'red')
    ax2.set_ylabel('fresh_water_flux',color = 'blue')
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')
    r = np.loadtxt(f'Data\{cell}\{cell}_r_daily2.txt')
    fresh_water_flux = np.loadtxt(f'Data\{cell}\{cell}_MB_daily2.txt')
    fresh_water_flux_brut = [] #Here we count only the fresh_water_flux period
    fresh_water_flux_net = [] # Here we count the net fresh_water_flux of the year (fresh_water_flux - freezing)
    for year in range(9):
        fresh_water_flux_brut_current = 0
        fresh_water_flux_net_current = 0
        for day in range(len(fresh_water_flux[year])):
            if not np.isnan(fresh_water_flux[year,day]):
                if fresh_water_flux[year,day] > 0 :
                    fresh_water_flux_brut_current += fresh_water_flux[year,day]
                fresh_water_flux_net_current += fresh_water_flux[year,day]
        fresh_water_flux_brut.append(fresh_water_flux_brut_current)
        fresh_water_flux_net.append(fresh_water_flux_net_current)
    ax1.plot(np.arange(2011,2020),r,color = 'red')
    #ax2.plot(np.arange(2011,2020),fresh_water_flux_brut,color = 'blue',linestyle = 'dashed')
    ax2.plot(np.arange(2011,2020),fresh_water_flux_net,color = 'blue')
    
    fig.suptitle(cell)
    ax1.grid()
    plt.show()
    plt.close()

#Mean fresh_water_flux
for cell in Cell_name:
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('r',color = 'red')
    ax2.set_ylabel('fresh_water_flux',color = 'blue')
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')
    r = np.loadtxt(f'Data\{cell}\{cell}_r_daily2.txt')
    fresh_water_flux = np.loadtxt(f'Data\{cell}\{cell}_MB_daily2.txt')
    fresh_water_flux_net = np.nanmean(fresh_water_flux,axis = 1)
    ax1.plot(np.arange(2011,2020),r,color = 'red')
    #ax2.plot(np.arange(2011,2020),fresh_water_flux_brut,color = 'blue',linestyle = 'dashed')
    ax2.plot(np.arange(2011,2020),fresh_water_flux_net,color = 'blue')
    
    fig.suptitle(cell)
    ax1.grid()
    plt.show()
    plt.close()