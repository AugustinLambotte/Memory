import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import fft, stats
import matplotlib as mpl

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
    kin = np.loadtxt(f'Data\{cell}\{cell}_index_daily.txt')
    fw_flux = np.loadtxt(f'Data\{cell}\{cell}_fw_daily.txt')
    siv = np.loadtxt(f'Data\{cell}\{cell}_siv_daily.txt')
    transp = np.loadtxt(f'Data/{cell}/{cell}_transp_daily.txt')

    #Convert data in good shape for interanual analysis
    kin_all = kin.flatten()[~np.isnan(kin.flatten())]
    fw_flux_all = fw_flux.flatten()[~np.isnan(fw_flux.flatten())]
    siv_all = siv.flatten()[~np.isnan(siv.flatten())]
    transp_all = transp.flatten()[~np.isnan(transp.flatten())]

    N = 365 #period in days of the running mean
    kin_run_mean = np.convolve(kin_all,np.ones(N)/N, mode = 'valid')
    fw_flux_run_mean = np.convolve(fw_flux_all,np.ones(N)/N, mode = 'valid')
    siv_run_mean = np.convolve(siv_all,np.ones(N)/N, mode = 'valid')
    transp_run_mean = np.convolve(transp_all,np.ones(N)/N, mode = 'valid')

    if N > 1: #If we do a running mean we lost the N first days of the record
        kin_run_mean = np.append(kin_run_mean,np.ones(N-1))
        fw_flux_run_mean = np.append(fw_flux_run_mean,np.ones(N-1))
        siv_run_mean = np.append(siv_run_mean,np.ones(N-1))
        transp_run_mean = np.append(transp_run_mean,np.ones(N-1))
        mask = np.concatenate((np.zeros(len(fw_flux_run_mean) - N + 1),(np.ones(N-1))))
        kin_run_mean = ma.masked_array(kin_run_mean, mask = mask)
        fw_flux_run_mean = ma.masked_array(fw_flux_run_mean, mask = mask)
        siv_run_mean = ma.masked_array(siv_run_mean, mask = mask)
        transp_run_mean = ma.masked_array(transp_run_mean, mask = mask)

    kin_all = kin_run_mean
    fw_flux_all = fw_flux_run_mean
    siv_all = siv_run_mean
    transp_all = transp_run_mean

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
    ax1.text(0.95, 0.07, f"r = {round(r_fw,2)}, p = {round(p_fw,3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right')
    ax2.set_ylim(-1,5)
    ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.arange(2011*365,2011*365+len(kin_all))/365,np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey')
    ax1.set_xlim(2011,2020)
    ax1.plot(np.arange(2011*365,2011*365+len(kin_all))/365,kin_all,color = 'red')
    ax2.plot(np.arange(2011*365,2011*365+len(fw_flux_all))/365,fw_flux_all,color = 'blue')
    if N > 1:
        fig.suptitle(f'Mean kinetic energy and fresh water flux over {cell} \n running average over {N} days', fontsize = 20)
    else:
        fig.suptitle(f'{cell} ', fontsize = 20)
    ax1.grid()
    plt.savefig(f'Plots/Correlation/4_cell/{cell}/{cell}_ALL_index_fw.png')
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
    ax1.text(0.95, 0.07, f"r = {round(r_siv,2)}, p = {round(p_siv,3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right')
    ax2.set_ylim(0,80)
    ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.arange(2011*365,2011*365+len(kin_all))/365,np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey')
    ax1.set_xlim(2011,2020)
    ax1.plot(np.arange(2011*365,2011*365+len(kin_all))/365,kin_all,color = 'red')
    ax2.plot(np.arange(2011*365,2011*365+len(fw_flux_all))/365,siv_all,color = 'blue')
    if N > 1:
        fig.suptitle(f'Mean kinetic energy and SIV over {cell} \n running average over {N} days', fontsize = 20)
    else:
        fig.suptitle(f'{cell} ', fontsize = 20)
    ax1.grid()
    plt.savefig(f'Plots/Correlation/4_cell/{cell}/{cell}_ALL_index_siv.png')
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
    ax1.text(0.95, 0.07, f"r = {round(r_transp,2)}, p = {round(p_transp,3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right')
    ax2.set_ylim(-1,5)
    ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.arange(2011*365,2011*365+len(kin_all))/365,np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey')
    ax1.set_xlim(2011,2020)
    ax1.plot(np.arange(2011*365,2011*365+len(kin_all))/365,kin_all,color = 'red')
    ax2.plot(np.arange(2011*365,2011*365+len(fw_flux_all))/365,transp_all,color = 'blue')
    if N > 1:
        fig.suptitle(f'Mean kinetic energy and Sea ice transport over {cell} \n running average over {N} days', fontsize = 20)
    else:
        fig.suptitle(f'{cell} ', fontsize = 20)
    ax1.grid()
    plt.savefig(f'Plots/Correlation/4_cell/{cell}/{cell}_ALL_index_transp.png')
    plt.close()
    

