import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import fft, stats
import matplotlib as mpl

"""
    Same as Cell_interannual_evol_bw.py but we compute the lag time with the maximum correlation

"""
mpl.rcParams['figure.figsize'] = (12,7)
mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['figure.titlesize'] = 20
mpl.rcParams['figure.labelsize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 15
mpl.rc('figure', titlesize = 50)
Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']


#Total fw_flux
for cell in Cell_name:
    #Extract data
    kin = np.loadtxt(f'Data/bw/{cell}\{cell}_index_daily.txt')
    fw_flux = np.loadtxt(f'Data/bw/{cell}\{cell}_fw_daily.txt')
    siv = np.loadtxt(f'Data/bw/{cell}\{cell}_siv_daily.txt')
    transp = np.loadtxt(f'Data/bw/{cell}/{cell}_transp_daily.txt')
    date = np.loadtxt('Data/bw/date.txt') #last columns are the day since 2011-01-01

    N = 20
    r_fw_lag = []
    p_fw_lag = []

    r_siv_lag = []
    p_siv_lag = []

    r_transp_lag = []
    p_transp_lag = []

    lag_time = [] # lag time value (nb of day)
    for lag in range(1,N):
        # Here, we compute r and p with a lag time between the kin_ene and the fw_flux from 0 days to N days.
        lag_time.append(date[lag,-1])
        kin_all = kin[abs(lag):]
        fw_flux_all = fw_flux[:-abs(lag)]
        siv_all = siv[:-abs(lag)]
        transp_all = transp[:-abs(lag)]

        r_fw_lag.append(stats.pearsonr(fw_flux_all,kin_all)[0])
        p_fw_lag.append(stats.pearsonr(fw_flux_all,kin_all)[1])

        r_siv_lag.append(stats.pearsonr(siv_all,kin_all)[0])
        p_siv_lag.append(stats.pearsonr(siv_all,kin_all)[1])

        r_transp_lag.append(stats.pearsonr(transp_all,kin_all)[0])
        p_transp_lag.append(stats.pearsonr(transp_all,kin_all)[1])
    
    plt.plot(range(len(r_fw_lag)),r_fw_lag)
    plt.plot(range(len(p_fw_lag)),p_fw_lag)
    plt.title('fw')
    plt.grid()
    plt.show()

    plt.plot(range(len(r_siv_lag)),r_siv_lag)
    plt.plot(range(len(p_siv_lag)),p_siv_lag)
    plt.title('siv')
    plt.grid()
    plt.show()

    plt.plot(range(len(r_transp_lag)),r_transp_lag)
    plt.plot(range(len(p_transp_lag)),p_transp_lag)
    plt.title('transp')
    plt.grid()
    plt.show()
    break

    r_fw = stats.pearsonr(kin_all,fw_flux_all)[0]
    p_fw = stats.pearsonr(kin_all,fw_flux_all)[1]

    
    r_siv = stats.pearsonr(kin_all,siv_all)[0]
    p_siv = stats.pearsonr(kin_all,siv_all)[1]

    r_transp = stats.pearsonr(kin_all,transp_all)[0]
    p_transp = stats.pearsonr(kin_all,transp_all)[1]

    ######## - FW flux - ###########

    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('J/kg',color = 'red')
    ax2.set_ylabel('km^3',color = 'blue')
    ax1.set_xticks(np.arange(2011,2021))
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')

    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax1.text(0.95, 0.95, f"r = {round(r_fw,2)}, p = {round(p_fw,3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right')
    #ax2.set_ylim(-1,5)
    #ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.linspace(2011,2020.541,len(kin_all)),np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey') #End date is 2020-07-16 which is 2020 + 54.1%
    #ax1.set_xlim(2011,2020)
    ax1.plot(np.linspace(2011,2020.541,len(kin_all)),kin_all,color = 'red') ################## ICI l'axe x va jusque 2020 mais les dernièrs données datent de 2020-07-16 donc faut CORRIGER
    ax2.plot(np.linspace(2011,2020.541,len(kin_all)),fw_flux_all,color = 'blue')
    
    fig.suptitle(f'Mean kinetic energy and net fresh water flux over {cell}', fontsize = 20)
    ax1.grid()
    plt.savefig(f'Plots/Correlation/4_cell/bw/{cell}/{cell}_ALL_index_fw.png')
    plt.close()

    ######## - SIV - ###########
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('J/kg',color = 'red')
    ax2.set_ylabel('km^3',color = 'blue')
    ax1.set_xticks(np.arange(2011,2021))
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')

    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax1.text(0.95, 0.95, f"r = {round(r_siv,2)}, p = {round(p_siv,3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right')
    #ax2.set_ylim(0,80)
    #ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.arange(len(kin_all)),np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey')
    #ax1.set_xlim(2011,2020)
    ax1.plot(np.arange(len(kin_all)),kin_all,color = 'red')
    ax2.plot(np.arange(len(fw_flux_all)),siv_all,color = 'blue')
    fig.suptitle(f'Mean kinetic energy and sea ice volume over {cell}', fontsize = 20)
    ax1.grid()
    plt.savefig(f'Plots/Correlation/4_cell/bw/{cell}/{cell}_ALL_index_siv.png')
    plt.close()

    ######## - Transport - ###########
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('J/kg',color = 'red')
    ax2.set_ylabel('km^3',color = 'blue')
    ax1.set_xticks(np.arange(2011,2021))
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')

    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax1.text(0.95, 0.95, f"r = {round(r_transp,2)}, p = {round(p_transp,3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right')
    #ax2.set_ylim(-1,5)
    #ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.arange(2011*365,2011*365+len(kin_all))/365,np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey')
    #ax1.set_xlim(2011,2020)
    ax1.plot(np.arange(2011*365,2011*365+len(kin_all))/365,kin_all,color = 'red')
    ax2.plot(np.arange(2011*365,2011*365+len(fw_flux_all))/365,transp_all,color = 'blue')
    fig.suptitle(f'Mean kinetic energy and net sea ice transport over {cell}', fontsize = 20)
    ax1.grid()
    plt.savefig(f'Plots/Correlation/4_cell/bw/{cell}/{cell}_ALL_index_transp.png')
    plt.close()
    

