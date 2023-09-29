import numpy as np
import matplotlib.pyplot as plt
from scipy import fft
Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']


#Total fw_flux
for cell in Cell_name:
    kin = np.loadtxt(f'Data\{cell}\{cell}_index_daily.txt')
    fw_flux = np.loadtxt(f'Data\{cell}\{cell}_MB_daily2.txt')
    r = np.loadtxt(f'Data\{cell}\{cell}_r_daily2.txt')
    p = np.loadtxt(f'Data\{cell}\{cell}_p_daily2.txt')
    for year in range(9):
        if any(np.isnan(kin[year])): #bissextile year
            kin_ft = fft.fftshift(fft.fft(kin[year],n = len(kin[year]) - 1))
            fw_flux_ft = fft.fftshift(fft.fft(fw_flux[year],n = len(fw_flux[year]) - 1))
        else:
            kin_ft = fft.fftshift(fft.fft(kin[year],n = len(kin[year])))
            fw_flux_ft = fft.fftshift(fft.fft(fw_flux[year],n = len(fw_flux[year])))
        kin_SP = np.abs(kin_ft)**2
        fw_flux_SP = np.abs(fw_flux_ft)**2
        plt.suptitle(f'Spectral power - {cell} - {year+2011}')
        plt.subplot(121)
        plt.title('Kinetic energy')
        plt.plot(np.arange(-len(kin_SP)/2,len(kin_SP)/2),kin_SP)
        plt.grid()
        plt.subplot(122)
        plt.title('Fresh water flux')
        plt.grid()
        plt.plot(np.arange(-len(fw_flux_SP)/2,len(fw_flux_SP)/2),fw_flux_SP)
        plt.savefig(f'Plots/Correlation/4_cell/{cell}/{year+2011}_Power_Spectrum.png')
        plt.clf()
        fig = plt.figure()
        ax1 = plt.axes()
        ax2 = plt.twinx()
        ax1.set_xlabel('month')
        ax1.set_ylabel('J/kg',color = 'red')
        ax2.set_ylabel('km^3',color = 'blue')
        ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
        ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')

        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        ax1.text(0.95, 0.07, f"r = {round(r[year],2)}, p = {round(p[year],3)}", fontsize=15, bbox=bbox,
            transform=ax1.transAxes, horizontalalignment='right')
        
        ax1.plot(range(len(kin[year,:])),kin[year,:],color = 'red')
        ax2.plot(range(len(fw_flux[year,:])),fw_flux[year,:],color = 'blue')
        
        fig.suptitle(f'{cell} - {year+2011}')
        ax1.grid()
        #plt.show()
        plt.savefig(f'Plots/Correlation/4_cell/{cell}/{year+2011}_index_and_MB.png')
        plt.close()