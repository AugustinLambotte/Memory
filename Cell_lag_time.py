import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, stats
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = (12,7)
mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['figure.titlesize'] = 40
mpl.rcParams['figure.labelsize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 15
Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']
cell_color = ['blue','orange','green','red']

#Total fw_flux
for cell,color in zip(Cell_name,cell_color):
    kin = np.loadtxt(f'Data\{cell}\{cell}_index_daily.txt')
    fw_flux = np.loadtxt(f'Data\{cell}\{cell}_MB_daily2.txt')
    
    
    #Convert data in good shape for interanual analysis
    kin_all = kin.flatten()[~np.isnan(kin.flatten())]
    fw_flux_all = fw_flux.flatten()[~np.isnan(fw_flux.flatten())]
    r,p = stats.pearsonr(fw_flux_all,kin_all)
    N = 2000
    r_lag = []
    p_lag = []
    for lag in range(1,N):
        # Here, we compute r and p with a lag time between the kin_ene and the fw_flux from 0 days to N days.

        kin = kin_all[abs(lag):]
        fw_flux = fw_flux_all[:-abs(lag)]
        r_lag.append(stats.pearsonr(fw_flux,kin)[0])
        p_lag.append(stats.pearsonr(fw_flux,kin)[1])
    plt.plot(np.arange(1,N)/365,r_lag,label = 'r',color = color)
    plt.grid()
    plt.ylim(-0.4,0.4)
    plt.title(f'{cell}',fontsize = 20)
    plt.yticks(color = color)
    plt.xticks([0,1,2,3,4,5])
    plt.ylabel('r', color = color)
    plt.xlabel('lag time (year)')
    plt.twinx()
    plt.ylim(0,5)
    plt.yticks([0,0.5,1],fontsize = 15, color = 'grey')
    plt.plot(np.arange(1,N)/365,p_lag, color = 'grey', label = 'p-value')
    plt.ylabel('p-value',y = 0.1,color = 'grey')
    #plt.legend()
    plt.savefig(f'Plots\Correlation/4_cell\{cell}/{cell}_lag_time.png')
    plt.clf()
        

