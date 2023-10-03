import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import fft
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
    fw_flux = np.loadtxt(f'Data\{cell}\{cell}_MB_daily2.txt')
    r = np.loadtxt(f'Data\{cell}\{cell}_r_daily2.txt')
    p = np.loadtxt(f'Data\{cell}\{cell}_p_daily2.txt')

    #Convert data in good shape for interanual analysis
    kin_all = kin.flatten()[~np.isnan(kin.flatten())]
    fw_flux_all = fw_flux.flatten()[~np.isnan(fw_flux.flatten())]
    r_all = r.flatten()[~np.isnan(r.flatten())]
    p_all = p.flatten()[~np.isnan(p.flatten())]

    N = 30 #period in days of the running mean
    kin_run_mean = np.convolve(kin_all,np.ones(N)/N, mode = 'valid')
    fw_flux_run_mean = np.convolve(fw_flux_all,np.ones(N)/N, mode = 'valid')
    if N > 1: #If we do a running mean we lost the N first days of the record
        kin_run_mean = np.insert(kin_run_mean,0,np.ones(N-1))
        fw_flux_run_mean = np.insert(fw_flux_run_mean,0,np.ones(N-1))
        mask = np.concatenate(((np.ones(N-1)),np.zeros(len(fw_flux_run_mean) - N + 1)))
        kin_run_mean = ma.masked_array(kin_run_mean, mask = mask)
        fw_flux_run_mean = ma.masked_array(fw_flux_run_mean, mask = mask)

    kin_all = kin_run_mean
    fw_flux_all = fw_flux_run_mean
    #Spectral analysis
    kin_ft = fft.fftshift(fft.fft(kin_all,n = len(kin_all)))
    fw_flux_ft = fft.fftshift(fft.fft(fw_flux_all,n = len(fw_flux_all)))
    kin_SP = np.abs(kin_ft)**2
    fw_flux_SP = np.abs(fw_flux_ft)**2
    if N > 1:
        plt.suptitle(f'Spectral power - {cell} - averaged over {N} days')
    else:
        plt.suptitle(f'Spectral power - {cell} ')
    plt.subplot(121)
    plt.title('Kinetic energy')
    plt.plot(np.arange(-len(kin_SP)/2,len(kin_SP)/2),kin_SP)
    plt.grid()
    plt.subplot(122)
    plt.title('Fresh water flux')
    plt.grid()
    plt.plot(np.arange(-len(fw_flux_SP)/2,len(fw_flux_SP)/2),fw_flux_SP)
    plt.savefig(f'Plots/Correlation/4_cell/{cell}/{cell}_ALL_year_Power_Spectrum.png')
    plt.clf()
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('J/kg',color = 'red')
    ax2.set_ylabel('km^3',color = 'blue')
    ax1.set_xticks(np.arange(2011,2021))
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')

    """ bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax1.text(0.95, 0.07, f"r = {round(r[year],2)}, p = {round(p[year],3)}", fontsize=15, bbox=bbox,
        transform=ax1.transAxes, horizontalalignment='right') """
    ax2.set_ylim(-4,15)
    ax1.set_ylim(0.004,0.0225)
    ax2.plot(np.arange(2011*365,2011*365+len(kin_all))/365,np.zeros(np.shape(kin_all)),linestyle = 'dashed',color = 'grey')
    ax1.set_xlim(2011,2020)
    ax1.plot(np.arange(2011*365,2011*365+len(kin_all))/365,kin_all,color = 'red')
    ax2.plot(np.arange(2011*365,2011*365+len(fw_flux_all))/365,fw_flux_all,color = 'blue')
    if N > 1:
        fig.suptitle(f'{cell} - running average over {N} days', fontsize = 20)
    else:
        fig.suptitle(f'{cell} ', fontsize = 20)
    ax1.grid()
    #plt.show()
    plt.savefig(f'Plots/Correlation/4_cell/{cell}/{cell}_ALL_year_index_and_MB.png')
    plt.close()
    

